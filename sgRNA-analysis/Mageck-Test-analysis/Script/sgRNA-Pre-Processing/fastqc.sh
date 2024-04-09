#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err

source /public3/soft/modules/module.sh
module load anaconda/3-Python-3.8.3-phonopy-phono3py
source activate py37

mkdir -p  ../result/mageck_count/fastqc
cd ../result/mageck_count/fastqc/

fastqc \
/public3/home/scg8403/Data/sgRNA/Data/PC_L1_1.fq.gz \
/public3/home/scg8403/Data/sgRNA/Data/PC_L1_2.fq.gz \
--outdir ./


fastqc \
/public3/home/scg8403/Data/sgRNA/Data/PG-1_L1_1.fq.gz \
/public3/home/scg8403/Data/sgRNA/Data/PG-1_L1_2.fq.gz \
--outdir ./

fastqc \
/public3/home/scg8403/Data/sgRNA/Data/PG-2_L1_1.fq.gz \
/public3/home/scg8403/Data/sgRNA/Data/PG-2_L1_2.fq.gz \
--outdir ./

fastqc \
/public3/home/scg8403/Data/sgRNA/Data/PG-3_L1_1.fq.gz \
/public3/home/scg8403/Data/sgRNA/Data/PG-3_L1_2.fq.gz \
--outdir ./

