#!/bin/bash

################################################################################
##################### GENERATE AMPLICONS WITH USEARCH ##########################
################################################################################

conda activate usearch-12.0-beta

inputdir="/home/storage/finished_projects/EvaG/ssUMI_amplicons_repeated_with_Illumina/raw_seq"
outputdir="/home/storage/finished_projects/EvaG/ssUMI_amplicons_repeated_with_Illumina/outidr"

# samples list (seperately for each primer)
ls -l | awk 'NR%2==1 {print $9}' > samples.list
sed -E -i 's/_L00[0-9]_R[0-9]_00[0-9]\.fastq//g' samples.list
nano samples.list

# Merge
while read sample; do
usearch -fastq_mergepairs $inputdir/v1v3/${sample}_L001_R1_001.fastq \
        -fastqout $sample.merged.fq \
        -threads 60 \
        -fastq_maxdiffs 15 \
        -sample ${sample}
cat $sample.merged.fq >> mergedv1v3.fq
done < $inputdir/v1v3/samples.list

while read sample; do
usearch -fastq_mergepairs $inputdir/v3v4/${sample}_L001_R1_001.fastq \
        -fastqout $sample.merged.fq \
        -threads 60 \
        -fastq_maxdiffs 15 \
        -sample ${sample}
cat $sample.merged.fq >> mergedv3v4.fq
done < $inputdir/v3v4/samples.list

# Strip primers
# V3-V4 (341F is 17 nt, 805R is 21 nt)
usearch -fastx_truncate mergedv3v4.fq -stripleft 17 -stripright 21 -padlen 0 -trunclen 400 -fastqout strippedv3v4.fq  # (1.9% too short)
# V1-V3 (27F is 20 nt, 534R is 17 nt)
usearch -fastx_truncate mergedv1v3.fq -stripleft 20 -stripright 17 -padlen 0 -trunclen 430 -fastqout strippedv1v3.fq # (3.2% too short)

# Quality filter
usearch -fastq_filter strippedv3v4.fq -fastq_maxee 1.0 -fastaout filteredv3v4.fa -relabel Filt -threads 60 # (81,3% passed)
usearch -fastq_filter strippedv1v3.fq -fastq_maxee 1.0 -fastaout filteredv1v3.fa -relabel Filt -threads 60 # (69.7% passed)

# Uniques
usearch -fastx_uniques filteredv3v4.fa -sizeout -relabel Uniq -minuniquesize 2 -fastaout uniquesv3v4.fa -threads $60 # (72.0% uniques)
usearch -fastx_uniques filteredv1v3.fa -sizeout -relabel Uniq -minuniquesize 2 -fastaout uniquesv1v3.fa -threads $60 # (75.2% uniques)

# Denoise
usearch -unoise3 uniquesv3v4.fa -zotus zotusv3v4.fa -relabel ZOTU -threads 60 -minsize 8 # default minsize (57,4 % chimeras)
usearch -unoise3 uniquesv1v3.fa -zotus zotusv1v3.fa -relabel ZOTU -threads 60 -minsize 8 # default minsize (40,2% chimeras)

# ZOTU table
usearch -otutab mergedv3v4.fq -otus zotusv3v4.fa -otutabout zotutab_rawv3v4.tsv -biomout ZOTUsv3v4.biom -threads 60  # (89.3% mapped)
usearch -otutab mergedv1v3.fq -otus zotusv1v3.fa -otutabout zotutab_rawv1v3.tsv -biomout ZOTUsv1v3.biom -threads 60  # (89,6% mapped)

################################################################################
########################### ASSIGN TAXONOMY WITH QIIME #########################
################################################################################

conda activate qiime2-amplicon-2024.10

# Import
qiime tools import \
    --input-path zotusv3v4.fa \
    --output-path zotusv3v4.qza \
    --type 'FeatureData[Sequence]'

qiime tools import \
      --input-path zotusv1v3.fa \
      --output-path zotusv1v3.qza \
      --type 'FeatureData[Sequence]'

# Vsearch
qiime feature-classifier classify-consensus-vsearch \
    --i-query zotusv3v4.qza \
    --i-reference-reads /home/storage/DB/silva_qiime2/silva-138-99-seqs.qza \
    --i-reference-taxonomy /home/storage/DB/silva_qiime2/silva-138-99-tax.qza \
    --p-perc-identity 0.97 \
    --o-classification vsearch_zotuclassification_v3v4.qza \
    --o-search-results vsearch_zotusearch_v3v4.qza

qiime feature-classifier classify-consensus-vsearch \
      --i-query zotusv1v3.qza \
      --i-reference-reads /home/storage/DB/silva_qiime2/silva-138-99-seqs.qza \
      --i-reference-taxonomy /home/storage/DB/silva_qiime2/silva-138-99-tax.qza \
      --p-perc-identity 0.97 \
      --o-classification vsearch_zotuclassification_v1v3.qza \
      --o-search-results vsearch_zotusearch_v1v3.qza
