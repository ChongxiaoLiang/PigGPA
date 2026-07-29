#!/bin/bash
#SBATCH -J gxe_model_t12_1
#SBATCH -p SMP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -t 00:30:00
#SBATCH -o /public/share/likui/hanyu/testdata/In-silico-data/t12/1/gxe_model.out
#SBATCH -e /public/share/likui/hanyu/testdata/In-silico-data/t12/1/gxe_model.err

cd /public/share/likui/hanyu/testdata/In-silico-data/t12/1
python3 gxe_model.py
