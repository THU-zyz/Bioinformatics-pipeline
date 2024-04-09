#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --output=mageck_P.out
#SBATCH --error=mageck_P.err

source /public3/soft/modules/module.sh
module load anaconda/3-Python-3.8.3-phonopy-phono3py
source activate py37

cd ../result/mageck_count/count
# 从fastq文件中收集sgRNA read count信息。输出count文件可以直接在test子命令中使用。
mageck count -l /public3/home/scg8403/Data/sgRNA/Data/sgRNA_library_Metabolic.csv -n PC9_all --sample-label "Treat1,Treat2,Treat3,Ctrl" --fastq /public3/home/scg8403/Data/sgRNA/Data/PG-1_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PG-2_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PG-3_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PC_L1.fq --pdf-report

mageck count -l /public3/home/scg8403/Data/sgRNA/Data/sgRNA_library_Metabolic.csv -n PC9_1 --sample-label "Treat1,Ctrl" --fastq /public3/home/scg8403/Data/sgRNA/Data/PG-1_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PC_L1.fq --pdf-report

mageck count -l /public3/home/scg8403/Data/sgRNA/Data/sgRNA_library_Metabolic.csv -n PC9_2 --sample-label "Treat2,Ctrl" --fastq /public3/home/scg8403/Data/sgRNA/Data/PG-2_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PC_L1.fq --pdf-report

mageck count -l /public3/home/scg8403/Data/sgRNA/Data/sgRNA_library_Metabolic.csv -n PC9_3 --sample-label "Treat3,Ctrl" --fastq /public3/home/scg8403/Data/sgRNA/Data/PG-3_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PC_L1.fq --pdf-report

# mageck test 用于比较sgRNA与基因的差异。

cd ../test
mageck test -k ../count/PC9_all.count.txt -t Treat1,Treat2,Treat3 -c Ctrl -n PC9_all --sort-criteria pos

mageck test -k ../count/PC9_1.count.txt -t Treat1 -c Ctrl -n PC9_1 --sort-criteria pos

mageck test -k ../count/PC9_2.count.txt -t Treat2 -c Ctrl -n PC9_2 --sort-criteria pos

mageck test -k ../count/PC9_3.count.txt -t Treat3 -c Ctrl -n PC9_3 --sort-criteria pos
