#!/bin/bash

read1=$1
read2=$2
sample=$3

cd /opt

mkdir out/$sample/$sample

mkdir out/$sample/fastqc
mkdir out/$sample/fastqc/trimmed
fastqc $1 $2 -o out/$sample/fastqc

cutadapt \
-q 18 \
--minimum-length 75 \
-o out/$sample/trimmed1.fq.gz \
-p out/$sample/trimmed2.fq.gz \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
-B CAAGCAGAAGACGGCATACGAGAT \
-O 15 \
-n 3 \
$1 $2

fastqc out/$sample/trimmed1.fq.gz out/$sample/trimmed2.fq.gz -o out/$sample/fastqc/trimmed

bbmerge.sh in1=out/$sample/trimmed1.fq.gz in2=out/$sample/trimmed1.fq.gz out=out/$sample/merged.fq.gz outu1=out/$sample/1_un.fq.gz outu2=out/$sample/2_un.fq.gz

rm out/$sample/trimmed*

bwa mem -M -R "@RG\tID:sample_name;cell;lane\tSM:sample_name\tPL:illumina\tLB:lib1" hologenome.fna out/$sample/1_un.fq.gz out/$sample/2_un.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > out/$sample/merged_un.bam

rm out/$sample/1_un.fq.gz
rm out/$sample/2_un.fq.gz

bwa mem -M -R "@RG\tID:sample_name;cell;lane\tSM:sample_name\tPL:illumina\tLB:lib1" hologenome.fna out/$sample/merged.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > out/$sample/merged.bam

rm out/$sample/merged.fq.gz

java -jar $PICARD MergeSamFiles I=out/$sample/merged.bam I=out/$sample/merged_un.bam SO=coordinate O=out/$sample/sorted_merged.bam

rm out/$sample/merged.bam
rm out/$sample/merged_un.bam

java -jar $PICARD MarkDuplicates \
REMOVE_DUPLICATES=true \
I=out/$sample/sorted_merged.bam \
O=out/$sample/dedup.bam \
M=out/$sample/mark_duplicates_report.txt \
VALIDATION_STRINGENCY=SILENT

rm out/$sample/sorted_merged.bam

samtools index out/$sample/dedup.bam

java -jar $GATK -T RealignerTargetCreator \
-R hologenome.fna \
-I out/$sample/dedup.bam \
-o out/$sample/dedup.intervals

java -jar $GATK \
-T IndelRealigner \
-R hologenome.fna \
-I out/$sample/dedup.bam \
-targetIntervals out/$sample/dedup.intervals \
-o out/$sample/dedup_indel.bam

rm out/$sample/dedup.bam*
rm out/$sample/dedup.intervals

samtools index out/$sample/dedup_indel.bam

samtools idxstats out/$sample/dedup_indel.bam > out/$sample/${sample}_duplicate_marked_idxstats.txt
samtools mpileup out/$sample/dedup_indel.bam -f hologenome.fna > out/$sample/${sample}_duplicate_marked_mpileup.txt

#Number of reads mapping to simulans and mel
grep -v "imulans" out/$sample/${sample}_duplicate_marked_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > output/$sample/num_mel.txt
grep "imulans" out/$sample/${sample}_duplicate_marked_idxstats.txt | awk -F '\t' '{sum+=$3;} END {print sum;}' > output/$sample/num_sim.txt

#Need to name these based on the sample
samtools view out/$sample/dedup_indel.bam | grep -v -P "imulans\t" > out/$sample/processed_mel_reads.bam
samtools view out/$sample/dedup_indel.bam | grep -P "imulans\t" > processed_sim_reads.bam

mv out/$sample/dedup_indel.bam  out/$sample/all_processed_reads.bam

python3 mpileup2sync.py --mpileup out/$sample/${sample}_duplicate_marked_mpileup.txt --ref dmel-all-chromosome-r6.12.fasta > out/$sample/output.sync
