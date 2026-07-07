# illumina_vs_nanopore_ssumi
The code in this repository is used to compare 16S rRNA gene amplicon sequencing reads with unique molecular identifiers on Oxford Nanopore, amplified with primer pair V1-V9 (ssUMI), to reads on Illumina, amplified with primer pairs V1-V3 and V3-V4.
The original ssUMI protocol was developed by Lin et al. (2024): https://doi.org/10.1093/pnasnexus/pgae411

Content:
- ssUMI_pipeline_and_amplicons.sh (merges UMI-tagged reads into consensus sequences, creates a ZOTU table, and assigns taxonomy to Nanopore reads)
- Illumina_amplicons (creates a ZOTU table and assigns taxonomy to Illumina reads)
- all_analysis.R (imports data into R; compares reads in two mock communities and five biological samples across different taxonomic levels)
