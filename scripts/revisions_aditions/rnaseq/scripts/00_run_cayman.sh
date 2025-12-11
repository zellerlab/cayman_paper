#!/bin/bash -e

# conda activate /g/typas/Personal_Folders/Nic/miniforge3/envs/cayman
cayman annotate_proteome ../data/v3 ../data/E_tayi_DSM26961.faa  --cutoffs ../data/v3/cutoffs.csv --output_file ../results/cayman_output/E_tayi_DSM26961.cayman --threads 8
cayman annotate_proteome ../data/v3 ../data/H_hathewayi_DSM13479.faa  --cutoffs ../data/v3/cutoffs.csv --output_file ../results/cayman_output/H_hathewayi_DSM13479.cayman --threads 8
