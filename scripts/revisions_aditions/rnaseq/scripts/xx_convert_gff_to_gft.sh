
#conda activate /g/scb2/zeller/karcher/mambaforge/envs/singularity
# Pull container
singularity pull docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0
# Start interactive session
singularity run --bind $(pwd)/../data  agat_1.4.2--pl5321hdfd78af_0.sif
# Run this then
for f in $(ls /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/data | grep gff$ | sed "s/\.gff//") ; do agat_convert_sp_gff2gtf.pl --gff /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/data/${f}.gff -o /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/data/${f}_from_agat.gtf ; done

