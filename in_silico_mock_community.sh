
# SRUN

# create taxon list with all species names and search for the species names in NCBI

datasets summary genome taxon \
  --inputfile taxon_list.txt \
  --assembly-level complete \
  --assembly-level chromosome \
  --assembly-source RefSeq \
  --as-json-lines | \
dataformat tsv genome \
  --fields \
accession,organism-name,organism-tax-id,organism-infraspecific-strain,organism-infraspecific-isolate,assminfo-biosample-strain,assminfo-biosample-collected-by,assminfo-biosample-isolation-source,ani-check-status,assminfo-sequencing-tech,assmstats-total-number-of-chromosomes,assmstats-number-of-contigs,assminfo-level,assminfo-atypicalis-atypical,assmstats-genome-coverage,assminfo-release-date,annotinfo-pipeline,annotinfo-release-date,assmstats-total-sequence-len,assmstats-contig-n50,assmstats-gc-percent,checkm-completeness,checkm-contamination,annotinfo-busco-complete,annotinfo-busco-lineage,annotinfo-featcount-gene-total,annotinfo-featcount-gene-protein-coding,type_material-display_text,assminfo-submitter,assminfo-biosample-accession,assminfo-bioproject \
> mock_community_genomes_metadata.tsv

# filter table in R
# import back into HPC

# download genomes

mkdir files

while read acc; do
  datasets download genome accession "$acc" \
    --include genome \
    --filename "files/${acc}.zip"
done < selected_accessions.tsv

# unzip genomes

mkdir temp_unzipped

for z in files/*.zip; do
name=$(basename "$z" .zip)
unzip -o "$z" -d "temp_unzipped/$name";
done

# move actual genomes out of subfolders

mkdir fasta_files

cp temp_unzipped/*/ncbi_dataset/data/*/*.fna fasta_files/

ls fasta_files/*.fna | wc -l

# use barrnap to extract rrna

mkdir rrnas

for file in fasta_files/*.fna; do
    name=$(basename "$file" .fna)
    barrnap \
      --kingdom bac \
      --outseq "rrnas/${name}.fa" \
      "$file" > "rrnas/${name}.gff"
done

# extract rrna coordinates

mkdir 16s_rrnas

for gff in rrnas/*.gff; do
name=$(basename "$gff" .gff)
awk -F'\t' '/Name=16S_rRNA/ {
printf "%s:%d-%d\n", $1, $4-1, $5
}' "$gff" > 16s_rrnas/${name}.txt
done

# use seqkit to search for 16S coordinates in files

for file in rrnas/*fa; do
name=$(basename "$file" .fa)
seqkit grep -n -r -f "16s_rrnas/${name}.txt" \
"$file" > "16s_rrnas/${name}.fa"
done

# count 16S per strain

grep -c ">" 16s_rrnas/*.fa | \
awk -F: '{print $2}' | \
sort -n | \
uniq -c

# check lenghth

cat 16s_rrnas/*.fa | \
seqkit fx2tab -n -l | \
awk '{print $2}' | \
sort -n | \
uniq -c

# use seqkit for n silico PCR V1-V9 (8F/1391R):

mkdir V1V9

for file in 16s_rrnas/*.fa; do
  name=$(basename "$file" .fa)
  seqkit amplicon "$file" \
  -F AGRGTTYGATYMTGGCTCAG \
  -R GACGGGCGGTGWGTRCA | \
  seqkit replace -p "^" -r "${name}|" \
  > "V1V9/${name}.fa"
done

# use seqkit for n silico PCR V1-V3 (27F/534R):

mkdir V1V3

for file in 16s_rrnas/*.fa; do
  name=$(basename "$file" .fa)
  seqkit amplicon "$file" \
  -F AGAGTTTGATCCTGGCTCAG \
  -R ATTACCGCGGCTGCTGG | \
  seqkit replace -p "^" -r "${name}|" \
  > "V1V3/${name}.fa"
done

# use seqkit for n silico PCR V3-V4 (341F/805R):

mkdir V3V4

for file in 16s_rrnas/*.fa; do
  name=$(basename "$file" .fa)
  seqkit amplicon "$file" \
  -F CCTACGGGNGGCWGCAG \
  -R GACTACHVGGGTATCTAATCC | \
  seqkit replace -p "^" -r "${name}|" \
  > "V3V4/${name}.fa"
done

# usearch loop

for region in V1V9 V1V3 V3V4; do

  echo "=== Processing $region ==="

  # Concatenate all per-genome fastas
  cat ${region}/*.fa > ${region}_all.fa

  # Add size=1 to each sequence (required by unoise3)
  awk '/^>/ {print $0 ";size=1"} !/^>/ {print}' ${region}_all.fa \
    > ${region}_all_sized.fa

  # Find unique sequences
  usearch -fastx_uniques ${region}_all_sized.fa \
    -sizeout \
    -relabel Uniq \
    -fastaout ${region}_uniques.fa \
    -threads 40

  # Denoise
  usearch -unoise3 ${region}_uniques.fa \
    -zotus ${region}_zotus.fa \
    -relabel ${region}_ZOTU \
    -minsize 1 \
    -threads 40

  # Map original sequences to ZOTUs
  usearch -otutab ${region}_all_sized.fa \
    -otus ${region}_zotus.fa \
    -otutabout ${region}_zotutab.tsv \
    -threads 40

done

# qiime Loop

for region in V1V9 V1V3 V3V4; do

  echo "=== Taxonomy for $region ==="

  # Import ZOTUs
  qiime tools import \
    --input-path ${region}_zotus.fa \
    --output-path ${region}_zotus.qza \
    --type 'FeatureData[Sequence]'

  # Assign taxonomy with vsearch
  qiime feature-classifier classify-consensus-vsearch \
    --i-query ${region}_zotus.qza \
    --i-reference-reads /home/storage/DB/silva_qiime2/silva-138-99-seqs.qza \
    --i-reference-taxonomy /home/storage/DB/silva_qiime2/silva-138-99-tax.qza \
    --p-perc-identity 0.97 \
    --o-classification ${region}_vsearchZOTUclassification.qza \
    --o-search-results ${region}_vsearchZOTUsearch.qza \
    --p-threads 40

  # Export taxonomy to tsv
  qiime tools export \
    --input-path ${region}_vsearchZOTUclassification.qza \
    --output-path ${region}_taxonomy/

done















































--------------------------------------------------------------------------------
# join and allign

cat V1V9/*.fa > V1V9_all.fa

mafft --auto --thread 40 V1V9_all.fa > V1V9_all_aln.fa

cat V1V3/*.fa > V1V3_all.fa

mafft --auto --thread 40 V1V3_all.fa > V1V3_all_aln.fa

cat V3V4/*.fa > V3V4_all.fa

mafft --auto --thread 40 V3V4_all.fa > V3V4_all_aln.fa

# use seqkit to tabulate

seqkit fx2tab V1V9_all.fa > V1V9_all_tab.txt

seqkit fx2tab V1V3_all.fa > V1V3_all_tab.txt

seqkit fx2tab V3V4_all.fa > V3V4_all_tab.txt

seqkit fx2tab V1V9_all_aln.fa > V1V9_all_aln_tab.txt

seqkit fx2tab V1V3_all_aln.fa > V1V3_all_aln_tab.txt

seqkit fx2tab V3V4_all_aln.fa > V3V4_all_aln_tab.txt
