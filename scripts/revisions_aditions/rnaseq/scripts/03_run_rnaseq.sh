#!/bin/bash -e

# To test, run
# nextflow run nf-core/rnaseq -profile test,singularity --outdir /scratch/karcher/an_outdir_lul -c ~/nextflow.confi

# Some notes
# nf-core/rnaseq is not designed for prokaryotic data, so we need some worksaround
# 0. Fix gtf: cat E_tayi_DSM26961.gtf |  grep -P '\btranscript_id\s+"[^"]+"'> E_tayi_DSM26961_no_transcript_id.gtf
# 1. generate an explicit transcript_fasta as such: gffread -w E_tayi_DSM26961_TRANSCRIPT.fna -g genomes/E_tayi_DSM26961.fna E_tayi_DSM26961_no_transcript_id.gtf
# 2. add --transcript_fasta E_tayi_DSM26961_TRANSCRIPT.fna and --extra_star_align_args "--sjdbGTFfeatureExon CDS" to the pipeline

mkdir -p ../results/nfcore_rnaseq/; conda activate /g/scb2/zeller/karcher/mambaforge/envs/nextflow; cd ../data/rnaseq_reads_proper_names_softlinks/;  nextflow run nf-core/rnaseq -w /scratch/karcher/c -with-timeline timeline_c -profile singularity --outdir ../../results/nfcore_rnaseq/H_hathewayi_DSM13479 -c ~/nextflow.confi --input ../../scripts/Hhathewayi_samplesheet.csv --fasta ../genomes/H_hathewayi_DSM13479.fna --gtf ../../data/H_hathewayi_DSM13479_no_transcript_id.gtf --transcript_fasta ../../data/H_hathewayi_DSM13479_TRANSCRIPT.fna  --save_align_intermeds --extra_star_align_args "--sjdbGTFfeatureExon CDS"; cd ../../scripts

mkdir -p ../results/nfcore_rnaseq/; conda activate /g/scb2/zeller/karcher/mambaforge/envs/nextflow; cd ../data/rnaseq_reads_proper_names_softlinks/;  nextflow run nf-core/rnaseq -w /scratch/karcher/d -with-timeline timeline_d -profile singularity --outdir ../../results/nfcore_rnaseq/E_tayi_DSM26961 -c ~/nextflow.confi --input ../../scripts/Etayi_samplesheet.csv --fasta ../genomes/E_tayi_DSM26961.fna --gtf ../../data/E_tayi_DSM26961_no_transcript_id.gtf --transcript_fasta ../../data/E_tayi_DSM26961_TRANSCRIPT.fna  --save_align_intermeds --extra_star_align_args "--sjdbGTFfeatureExon CDS"  --save_align_intermeds; cd ../../scripts
