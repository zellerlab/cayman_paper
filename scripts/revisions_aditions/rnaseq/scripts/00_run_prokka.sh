#conda activate prokka_env
export PERL5LIB="/g/typas/Personal_Folders/Nic/miniforge3/envs/prokka_env/lib/perl5/site_perl/"
mkdir -p ../results/prokka
#for f in $(ls ../results/assemblies_softlinks/ | sed "s/\.fasta//")
#for f in $(ls ../results/assemblies_softlinks_selected/ | sed "s/\.fasta//")
for f in $(ls ../data/genomes/ | grep fna$ | sed "s/\.fna//")
do
        echo -e "prokka --cpus 4 --force --outdir ../results/prokka/${f} --prefix ${f} ../data/genomes/${f}.fna"
done > iv.sh
