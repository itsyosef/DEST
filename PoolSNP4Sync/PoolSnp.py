import sys
from collections import defaultdict as d
import re
from optparse import OptionParser, OptionGroup
import math
import gzip
from datetime import datetime

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="""python %prog \
      --sync data.sync \
      --min-cov 10 \
      --max-cov data.cov \
      --min-count 10 \
      --min-freq 0.01 \
      --mis-frac 0.1 \
      --names sample1,sample2 \
      > output.vcf"""
parser = OptionParser(usage=usage)
helptext="""

H E L P :
_________
"""
group=OptionGroup(parser,helptext)
#########################################################   parameters   #########################################################################

parser.add_option("--sync", dest="m", help="A sync file")
parser.add_option("--min-cov", dest="minc", help="The minimum coverage threshold: e.g. 10",default=10)
parser.add_option("--max-cov", dest="max", help="An input file with precomputed coverage thresholds")
parser.add_option("--min-count", dest="mint", help="The minimum number of counts of the alternative allele across all samples pooled",default=3)
parser.add_option("--min-freq", dest="minf", help="The minimum Frequency of the alternative allele across all samples pooled",default=0.01)
parser.add_option("--miss-frac", dest="mis", help="The minimum Frequency of the alternative allele across all samples pooled",default=0.1)
parser.add_option("--names", dest="n", help="a comma separted list of thenames of all samples in the sync file")

parser.add_option_group(group)
(options, args) = parser.parse_args()


################################### functions ######################################


def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''

  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

def sync2string(x):
    ''' convert sync format to string of nucleotides  where x is a string in sync format '''
    string=""
    if x==".:.:.:.:.:." or x=="0:0:0:0:0:0":
        return "na"
    alleles=["A","T","C","G"]
    ah=dict(zip(alleles,map(int,x.split(":")[:4])))
    for k,v in ah.items():
        string+=v*k
    return string
covh=d(lambda:d(int))

def keywithmaxvalue(x):
    ''' This function resturns the key for the maximum value in a dictionary'''
    newhash=d(list)
    for k,v in x.items():
        newhash[v].append(k)
    return newhash[max(newhash.keys())]

################################## parameters ########################################

data=options.m
minimumcov=int(options.minc)
minimumcount=int(options.mint)
minimumfreq=float(options.minf)
missfrac=float(options.mis)

############################ get MAX coverage threshold  #############################
maximumcov=d(list)
for l in open(options.max,"r"):
    if l.startswith("#") or l.startswith("calculating"):
        continue
    k,v=l.split("\t")
    maximumcov[k]=[int(x) for x in v.split(",")]
#print maximumcov
############################ parse sync ###########################################

# parse sync and store alternative alleles:

# Returns a datetime object containing the local date and time
dateTimeObj = datetime.now()

print("##fileDate="+str(dateTimeObj.day)+"/"+str(dateTimeObj.month)+"/"+str(dateTimeObj.year))
print("##Source=PoolSnp")
print("##Parameters=<ID=MinCov,Number="+options.minc+",Type=Integer,Description=\"Minimum coverage per sample\">")
print("##Parameters=<ID=MinCount,Number="+options.mint+",Type=Integer,Description=\"Minimum alternative allele count across all samples pooled\">")
print("##Parameters=<ID=MinFreq,Number="+options.minf+",Type=Integer,Description=\"Minimum alternative allele frequency across all samples pooled\">")
print("##Parameters=<ID=MaximumMissingFraction,Number="+options.mis+",Type=Float,Description=\"Maximum fraction of samples allowed that are not fullfilling all parameters\">")
print("""##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average per-sample depth of bases\">
##INFO=<ID=NC,Number=1,Type=Integer,Description=\"Number of samples not called\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Reference Counts\">
##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Alternative Counts\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##FORMAT=<ID=FREQ,Number=1,Type=FLoat,Description=\"Variant allele frequency\">""")
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(options.n.split(",")))

for l in load_data(data):
    a=l.rstrip().split()

    CHR,POS,REF = a[:3]

    ## only keep chromosomal arms with maximum coverage threshold
    if CHR not in maximumcov:
        #print CHR
        continue

    libraries=a[3:]
    # loop through libraries
    totalalleles=d(int)
    alleles=d(lambda:d(int))
    covtest=d(int)

    for j in range(len(libraries)):
        alleles[j]
        nuc = sync2string(libraries[j])

        # test if seq-string is empty
        if nuc=="na":
            covtest[j]=1
            continue

        # ignore if coverage is below or above thresholds after filtering for 1) InDels and 2) base-quality
        if len(nuc)<minimumcov or len(nuc)>maximumcov[CHR][j]:
            covtest[j]=1
            continue
        else:
            covtest[j]=0

        # read all alleles
        for i in range(len(nuc)):

            # count alternative nucleotides
            totalalleles[nuc[i].upper()]+=1
            alleles[j][nuc[i].upper()]+=1

    ## test if SNPs pass minimum count / minimum frequency threshold:
    for allele,counts in totalalleles.items():
        if counts<minimumcount or counts/float(sum(totalalleles.values()))<minimumfreq:
            del totalalleles[allele]

    ## test if site is polymorphic
    if len(totalalleles)<2:
        #print CHR,POS,"non-poly",totalalleles
        continue

    ## create output for VCF
    ADP=sum(totalalleles.values())/len(libraries)
    ALT=[]
    ## set alternative allele order:
    for i in ["A","T","C","G"]:
        if i==REF:
            continue
        if i not in totalalleles:
            continue
        ALT.append(i)


    ## set ADP,NC,GT,AD and DP
    ADP=sum(totalalleles.values())/len(libraries)
    samplelist=[]
    co=0
    miss=0

    for j in range(len(libraries)):
        ## make empty entry if no allele counts for sample
        if j not in alleles:
            samplelist.append("./.:.:.:.:.")
            miss+=1
            continue

        ## make empty entry if sample not fullfilling min/max coverage threshold
        if covtest[j]==1:
            samplelist.append("./.:.:.:.:.")
            miss+=1
            continue

        alleleh = alleles[j]
        # remove alleles not counted in all samples
        for k,v in alleleh.items():
            if k != REF and k not in ALT:
                del alleleh[k]
        GT,AD,RD,FREQ,NC=[],[],0,[],0
        DP=sum(alleleh.values())

        ## test if mincoverage is still reached when removing alleles that do not fullfill criteria
        if DP<minimumcov:
            samplelist.append("./.:.:.:.:.")
            continue

        ## test if sample empty:
        if len(alleleh)==0:
            NC+=1
            samplelist.append("./.:.:.:.:.")
            miss+=1
            continue

        # test if population is fixed for REF allele
        if len(alleleh)==1 and REF in alleleh:
            samplelist.append("0/0:"+str(DP)+":0:"+str(DP)+":0.0")
            continue

        # test if population is fixed for ALT allele
        at=0
        if len(alleleh)==1:
            for i in range(len(ALT)):
                if ALT[i] in alleleh:
                    samplelist.append(str(i+1)+"/"+str(i+1)+":0:"+str(DP)+":"+str(DP)+":1.0")
                    at=1
                    continue
        if at==1:
            continue

        ## proceed if population not fixed
        ## set REF counts
        if REF in alleleh:
            GT.append(0)
        ## set ALT counts
        for i in range(len(ALT)):
            if ALT[i] in alleleh:
                GT.append(i+1)
                AD.append(alleleh[ALT[i]])
                RD=DP-sum(AD)
                FREQ.append(round(alleleh[ALT[i]]/float(sum(alleleh.values())),2))

        samplelist.append("/".join(map(str,GT))+":"+str(RD)+":"+",".join(map(str,AD))+":"+str(DP)+":"+",".join(map(str,FREQ)))

    ## test if missing fraction of samples smaller than threshold:
    if miss/float(len(libraries))>missfrac:
        #print CHR,POS,"missing fraction",miss/float(len(libraries))
        continue


    ## write output
    print(CHR+"\t"+POS+"\t.\t"+REF+"\t"+",".join(ALT)+"\t.\t.\tADP="+str(ADP)+";NC="+str(NC)+"\tGT:RD:AD:DP:FREQ\t"+"\t".join(samplelist))
