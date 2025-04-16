
#conda activate /g/scb2/zeller/karcher/mambaforge/envs/singularity
# Pull container
singularity pull docker://quay.io/biocontainers/agat:1.4.2--pl5321hdfd78af_0
# Start interactive session
singularity run --bind $(pwd)/../results/prokka/ agat_1.4.2--pl5321hdfd78af_0.sif
# Run this then
for f in $(ls /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/prokka/ | grep -E "C_|E_|H_|P_|B_"); do agat_convert_sp_gff2gtf.pl --gff /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/prokka/${f}/${f}.gff -o /g/scb/zeller/karcher/cayman_paper/scripts/revisions_aditions/rnaseq/results/prokka/${f}/${f}.gtf ; done
