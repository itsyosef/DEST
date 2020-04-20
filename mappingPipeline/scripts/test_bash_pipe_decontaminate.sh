#!/bin/bash

mkdir out

java -jar $PICARD CreateSequenceDictionary \
R=decontaminate.fasta \
O=decontaminate.dict

samtools faidx decontaminate.fasta

mkdir out/fastqc
mkdir out/fastqc/trimmed
fastqc $1 $2 -o fastqc

cutadapt \
-q 18 \
--minimum-length 75 \
-o out/sub_1_trim.fq.gz \
-p out/sub_2_trim.fq.gz \
-b ACACTCTTTCCCTACACGACGCTCTTCCGATC \
-B CAAGCAGAAGACGGCATACGAGAT \
-O 15 \
-n 3 \
$1 $2

fastqc out/sub_1_trim.fq.gz out/sub_2_trim.fq.gz -o out/fastqc/trimmed

bbmerge.sh in1=out/sub_1_trim.fq.gz in2=out/sub_2_trim.fq.gz out=out/merged.fq.gz outu1=out/sub_1_un.fq.gz outu2=out/sub_2_un.fq.gz

bwa mem -M -R "@RG\tID:sample_name;cell;lane\tSM:sample_name\tPL:illumina\tLB:lib1" decontaminate.fasta out/sub_1_un.fq.gz out/sub_2_un.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > out/merged_un.bam

bwa mem -M -R "@RG\tID:sample_name;cell;lane\tSM:sample_name\tPL:illumina\tLB:lib1" decontaminate.fasta out/merged.fq.gz | samtools view -Sbh -q 20 -F 0x100 - > out/merged.bam

java -jar $PICARD MergeSamFiles I=out/merged.bam I=out/merged_un.bam SO=coordinate O=out/sorted_merged.bam

java -jar $PICARD MarkDuplicates \
REMOVE_DUPLICATES=true \
I=out/sorted_merged.bam \
O=out/dedup.bam \
M=out/dedup_report.txt \
VALIDATION_STRINGENCY=SILENT

samtools index out/dedup.bam

java -jar $GATK -T RealignerTargetCreator \
-R decontaminate.fasta \
-I out/dedup.bam \
-o out/dedup.intervals

java -jar $GATK \
-T IndelRealigner \
-R decontaminate.fasta \
-I out/dedup.bam \
-targetIntervals out/dedup.intervals \
-o out/dedup_indel.bam

samtools idxstats out/dedup_indel.bam > out/dedup_idxstats.txt
samtools mpileup out/dedup_indel.bam -f decontaminate.fasta > out/dedup_mpileup.txt

#Need to name these based on the reference
samtools view dedup_indel.bam | grep -v -P "imulans\t" > processed_mel_reads.bam
samtools view dedup_indel.bam | grep -P "imulans\t" > processed_sim_reads.bam

python3 mpileup2sync.py --mpileup out/dedup_mpileup.txt --ref dmel-all-chromosome-r6.12.fasta > output.sync
