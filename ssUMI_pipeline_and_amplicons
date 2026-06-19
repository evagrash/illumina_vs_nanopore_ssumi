#!/bin/bash

################################################################################
################################ MERGE RAW FILES ###############################
################################################################################

# Locate raw files
cd home/storage/ssUMI_amplicons_NLZOH_2025-08-18

# Make output folder
mkdir /home/storage/finished_projects/EvaG/ssUMI_amplicons_NLZOH_2025-08-18_merged

# For each flowcell, merge all passed reads for each barcode
for folder in /home/storage/finished_projects/EvaG/ssUMI_amplicons_NLZOH_2025-08-18/ssUMI_amplicon_LIB2/20250818_1402_X2_FBC90119_1f64a87e/fastq_pass/*; do
  barcodename=$(basename "$folder")
  cat "$folder"/*.fastq.gz > /home/storage/finished_projects/EvaG/ssUMI_amplicons_NLZOH_2025-08-18_merged/"$barcodename"_L1.gz
done

# Unzip files (keep originals)
for file in *.gz; do
    gunzip -c "$file" > "${file%.gz}"
done

# Rename files using a list
while IFS=';' read oldname newname ; do
  mv "$oldname" "$newname" ;
done < rename_nanopore_samples.csv

################################################################################
########################## START ssUMI PIPELINE ################################
################################################################################

# https://github.com/ZielsLab/ssUMI
# Install ssUMI pipeline as described on GitHUb
# Here, the rapid mode will be used, as it does not require a GPU
# To use ssUMI in standard mode, replace ssumi_rapid with ssumi_std, add -n and -T values

# Input file is a signle file containing raw Nanopore data in fastq format
# Primers -F and -R anneal to 16S
# Adaptors -f and -r are what connects the UMI tags (ssUMI-8F-UMI/ssUMI-1391R-UMI) and universal UMI multipying primers (ssUMI-Universal-F/ssUMI-Universal-R)
# Medaka models are named after: {pore}_{device}_{caller variant}_{caller version}

srun -N 1 --mem=100GB --cpus-per-task=200 --pty /bin/bash

for file in *.fastq; do
    name=$(basename "$file")
    longread_umi ssumi_rapid \
      -d $file \
      -v 3 \
      -o "$name"_after_UMI \
      -s 200 \
      -e 200 \
      -E 0.1 \
      -m 1200 \
      -M 2000 \
      -f GTATCGTGTAGAGACTGCGTAGG \
      -F AGRGTTYGATYMTGGCTCAG \
      -r AGTGATCGAGTCAGTGCGAGTG \
      -R GACGGGCGGTGWGTRCA \
      -c 3 \
      -q r1041_e82_400bps_hac_g632 \ #medaka model
      -t 40

############################################################################################
#################################### IMPORT INTO USEARCH ###################################
############################################################################################

eval "$(conda shell.bash hook)"
source activate usearch12.0-beta

# Create output file
> all_samples.fasta

# Rename fasta headers from ssUMi format [>umi3733bins;ubs=4] to usearch format: [Sample.37333 37333 ubs=5]
for seq in */consensus_raconx3.fa; do

  folder=$(basename "$(dirname "$seq")")
  name=${folder%_after_UMI}

  awk -v s="$name" '
              /^>/ {
              match($0,/umi([0-9]+)/,a)
              number=a[1]
              match($0,/ubs=[0-9]+/,b)
              coverage=b[0]
              print ">"s"."number" "number" "coverage
              next
          }
          {print}
          ' "$seq" >> all_samples.fasta
done

################################################################################
############################### WRITE OUT UMI COUNTS ###########################
################################################################################

# Write counts into file
echo -e "Sample_Name\tubs_count" > umi_counts.tsv

# Format
grep "^>" all_samples.fasta | awk -F'[. ]' '{
    sample=$1; gsub(">", "", sample);
    count=$4; gsub("ubs=", "", count);
    print sample "\t" count
}' >> umi_counts.tsv

################################################################################
########################## MAKE ZOTUWITH USEARCH ##############################
################################################################################

# Find unique read counts and abundances
usearch -fastx_uniques all_samples.fasta -sizeout -relabel Uniq -fastaout uniques.fa -threads 40

# Denoise: predict biological sequences and filter chimeras
usearch -unoise3 uniques.fa -zotus zotus.fa -relabel ZOTU -minsize 8 -threads 40

# Make ZOTU table
usearch -otutab all_samples.fasta -otus zotus.fa -otutabout zotutab_raw.tsv -biomout ZOTUs.biom -threads 40

################################################################################
######################### PREDICT TAXONOMY WITH QIIME ##########################
################################################################################

eval "$(conda shell.bash hook)"
source activate qiime2-amplicon-2024.10

# Import
qiime tools import \
    --input-path zotus.fa \
    --output-path zotus.qza \
    --type 'FeatureData[Sequence]'

# Assign taxonomy with vsearch
qiime feature-classifier classify-consensus-vsearch \
        --i-query zotus.qza \
        --i-reference-reads /home/storage/DB/silva_qiime2/silva-138-99-seqs.qza \
        --i-reference-taxonomy /home/storage/DB/silva_qiime2/silva-138-99-tax.qza \
        --p-perc-identity 0.97 \
        --o-classification vsearchZOTUclassification.qza \
        --o-search-results vsearchZOTUsearch.qza
