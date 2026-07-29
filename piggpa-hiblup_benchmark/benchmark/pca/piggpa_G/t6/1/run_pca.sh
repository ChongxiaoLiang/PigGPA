#!/bin/bash
#SBATCH -J pca_t6_1
#SBATCH -p SMP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:30:00
#SBATCH -o /public/share/likui/hanyu/testdata/In-silico-data/t6/1/pca_analysis.out
#SBATCH -e /public/share/likui/hanyu/testdata/In-silico-data/t6/1/pca_analysis.err

cd /public/share/likui/hanyu/testdata/In-silico-data/t6/1
python3 pca_analysis.py
