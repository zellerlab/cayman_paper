#!/bin/bash -e

# To test, run
# nextflow run nf-core/rnaseq -profile test,singularity --outdir /scratch/karcher/an_outdir_lul -c ~/nextflow.confi

mkdir -p ../results/nfcore_rnaseq/; conda activate /g/scb2/zeller/karcher/mambaforge/envs/nextflow; cd ../data/rnaseq_reads_proper_names_softlinks/;  nextflow run nf-core/rnaseq -w /scratch/karcher/c -with-timeline timeline_c -profile singularity --outdir ../../results/nfcore_rnaseq/H_hathewayi_DSM13479 -c ~/nextflow.confi --input ../../scripts/Hhathewayi_samplesheet.csv --fasta ../genomes/H_hathewayi_DSM13479.fna --gtf ../../data/H_hathewayi_DSM13479_from_agat.gtf  --skip_gtf_transcript_filter  --skip_dupradar --save_align_intermeds; cd ../../scripts

mkdir -p ../results/nfcore_rnaseq/; conda activate /g/scb2/zeller/karcher/mambaforge/envs/nextflow; cd ../data/rnaseq_reads_proper_names_softlinks/;  nextflow run nf-core/rnaseq -w /scratch/karcher/d -with-timeline timeline_d -profile singularity --outdir ../../results/nfcore_rnaseq/E_tayi_DSM26961 -c ~/nextflow.confi --input ../../scripts/Etayi_samplesheet.csv --fasta ../genomes/E_tayi_DSM26961.fna --gtf ../../data/E_tayi_DSM26961_from_agat.gtf --skip_gtf_transcript_filter  --skip_dupradar --save_align_intermeds; cd ../../scripts
