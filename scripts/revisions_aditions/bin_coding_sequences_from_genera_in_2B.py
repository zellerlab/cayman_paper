import pandas as pd
import os
import numpy as np
import math
from Bio import SeqIO

outputfolderbase = '/g/scb/zeller/karcher/cayman_paper/data/'
os.makedirs(outputfolderbase, exist_ok = True)

data = pd.read_csv('../../data/almeida_cazy_annotations.tsv', sep = "\t")
goi = [
    "Eisenbergiella",
    "Hungatella",
    "Akkermansia",
    "Bacteroides",
    "Barnesiella",
    "Coprobacter",
    "Parabacteroides",
    "Paraprevotella"
]
goi = pd.DataFrame(goi, columns = ['Genus'])
# fix genus column
data = data[[True if (isinstance(x, str)  or not math.isnan(x)) else False for x in data['Genus']]]
data = data[[True if not x.startswith('NA') else False for x in data['Genus']]]
data['Genus'] = [x.split(' ')[1] for x in data['Genus']]
# filter for genera of interest
data = pd.merge(data, goi, on = 'Genus', how = 'inner')
# Loop over entries in the genomeID column and create softlinks by genus

with open('/g/scb2/zeller/SHARED/DATA/assembled_genomes/Almeida_2020_combined_set/find.prokka.whole.path', "r") as f:
    paths = f.readlines()
    paths = [x.strip() for x in paths]

genomes = [path.split('/')[-1].replace(".faa", "") for path in paths]
paths = pd.DataFrame([genomes, paths], index = ['genome', 'Path']).T

from collections import defaultdict
orfid_cazyfam_map = defaultdict(list)
#orfid_cazyfam_map = dict(zip(data['sequenceID'], data['cazy_family']))
for orfid, cazyfam in zip(data['sequenceID'], data['cazy_family']):
    orfid_cazyfam_map[orfid].append(cazyfam)

for genus in goi['Genus']:
    print(f"Working on {genus}")
    os.makedirs(os.path.join(outputfolderbase, 'almeida_cazy_info_by_genus/', genus), exist_ok = True)
    os.makedirs(os.path.join(outputfolderbase, 'almeida_coding_sequences_by_genus/', genus), exist_ok = True)
    data_genus = data[data['Genus'] == genus]
    # Write genus-specific cazy information for quinten
    data_genus.to_csv(os.path.join(outputfolderbase, 'almeida_cazy_info_by_genus/', genus, f'{genus}_cazy_annotations.tsv'), sep = "\t", index = False)
    data_genus = data_genus[['genome']]
    data_genus = data_genus.drop_duplicates()
    paths_genus = pd.merge(paths, data_genus, how = 'inner', on = "genome")
    # Remove all files in os.path.join(outputfolderbase, genus)
    dev_null = [os.remove(os.path.join(outputfolderbase, 'almeida_coding_sequences_by_genus/', genus, x)) for x in os.listdir(os.path.join(outputfolderbase, "almeida_coding_sequences_by_genus", genus))]
    for genome, path_base in zip(paths_genus['genome'], paths_genus['Path']):
        print(genome)
        path_base = path_base.replace('.faa', '')
        faa = list(SeqIO.parse(path_base + ".faa", "fasta"))
        # Very basic gff parsing... error prone but should be fine for this purpose
        with open(path_base + ".gff", "r") as f:
            gff = f.readlines()
            gff = [x for x in gff if x.startswith('GUT')]
            gff = [x.split('\t') for x in gff]
            gff = pd.DataFrame(gff)
            gff['ID'] = [x.split(';')[0].split('=')[1] for x in gff[8]]
            gff['seq_id_out'] = ["___".join([a,b,c,d,e,f]) for a,b,c,d,e,f in zip(gff[0], gff[3], gff[4], gff[6], gff['ID'], [";;;".join(orfid_cazyfam_map[id]) for id in gff['ID']])]
            gff = dict(zip(gff['ID'], gff['seq_id_out']))
            for record in faa:
                record.id = gff[record.id]
        with open(os.path.join(outputfolderbase, 'almeida_coding_sequences_by_genus/', genus, genome + '.faa'), "w") as f_out:
            SeqIO.write(faa, f_out, "fasta")


    #dev_null = [os.symlink(path, os.path.join(outputfolderbase, 'almeida_coding_sequences_by_genus/', genus, genome + '.faa')) for genome, path in zip(paths_genus['genome'], paths_genus['Path'])]
    # Also link gffs...
    #dev_null = [os.symlink(path.replace('.faa', '.gff'), os.path.join(outputfolderbase, 'almeida_coding_sequences_by_genus/', genus, genome + '.gff')) for genome, path in zip(paths_genus['genome'], paths_genus['Path'])]

    

