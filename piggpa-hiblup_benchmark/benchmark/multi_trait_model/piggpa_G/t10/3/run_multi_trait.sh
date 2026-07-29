#!/bin/bash
#SBATCH -J multi_trait_t10_3
#SBATCH -p SMP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:30:00
#SBATCH -o /public/share/likui/hanyu/testdata/In-silico-data/t10/3/multi_trait.out
#SBATCH -e /public/share/likui/hanyu/testdata/In-silico-data/t10/3/multi_trait.err

cd /public/share/likui/hanyu/testdata/In-silico-data/t10/3
python3 multi_trait_model.py
