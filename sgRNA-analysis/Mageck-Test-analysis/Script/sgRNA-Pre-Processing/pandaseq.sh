#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --output=pandaseq.out
#SBATCH --error=pandaseq.err

source /public3/soft/modules/module.sh
module load anaconda/3-Python-3.8.3-phonopy-phono3py
source activate py37

cd /public3/home/scg8403/Data/sgRNA/Data
# 合并双端测序fastq文件
pandaseq -F -f PG-1_L1_1.fq \
-r PG-1_L1_2.fq \
-w PG-1_L1.fq

pandaseq -F -f PG-2_L1_1.fq \
-r PG-2_L1_2.fq \
-w PG-2_L1.fq

pandaseq -F -f PG-3_L1_1.fq \
-r PG-3_L1_2.fq \
-w PG-3_L1.fq



