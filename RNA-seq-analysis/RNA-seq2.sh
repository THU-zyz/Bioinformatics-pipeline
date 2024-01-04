#!/bin/bash
#SBATCH -p amd_512
#SBATCH -n 1
#SBATCH -c 32

source /public3/soft/modules/module.sh
module load anaconda/3-Python-3.8.3-phonopy-phono3py
source activate py37


##### settings #####
# 如样本名为 8-cell_re1_R1.fastq.gz
# 样本名 如 "8-cell_re1"
sample=$1
# 样本所在目录 如 "/public3/home/scg8403/Data/人胚胎RNA-seq/Data"
path=$2
# 输出结果文件夹 如 "/public3/home/scg8403/Data/人胚胎RNA-seq/result"
output_path=$3
# mapping 的基因组目录路径  如 "/public3/home/scg8403/Data/ref_genome/index/STAR/genome"
genome_index_path=$4
# rRNA基因组目录路径 如 "/public3/home/scg8403/Data/ref_genome/index/rRNA/Homo_sapiens.rRNA"
rRNA_index=$5
# 基因组GTF文件所在位置 如 "/public3/home/scg8403/Data/ref_genome/index/GTF/Homo_sapiens.GRCh38.110.gtf"
GTF=$6
##### settings #####



# 1. fastqc 质控检测
#echo start fastqc ${sample} `date`
#cd ${output_path}
#mkdir -p fastqc/${sample}

#fastqc \
#${path}/${sample}_R1.fq.gz \
#${path}/${sample}_R2.fq.gz \
#--outdir ${output_path}/fastqc/${sample} --noextract >> ${output_path}/fastqc/${sample}/${sample}.log 2>&1

#echo finish fastqc ${sample} `date`

# 2. fastp 除去低质量片段
echo start fastp ${sample} `date`
cd ${output_path}
mkdir -p quality_control/${sample}

fastp -i ${path}/${sample}_R1.fq.gz \
-I ${path}/${sample}_R2.fq.gz \
-o ${output_path}/quality_control/${sample}/${sample}.clean.1.fastq.gz \
-O ${output_path}/quality_control/${sample}/${sample}.clean.2.fastq.gz \
--thread=32 -l 15 \
-j ${output_path}/quality_control/${sample}/${sample}.json \
-h ${output_path}/quality_control/${sample}/${sample}.html


echo finish $fastp ${sample} `date`

# 3. bowtie rRNA参考基因组以去除 rRNA部分
echo start remove rRNA ${sample} `date`
cd ${output_path}
mkdir -p remove_rRNA/fastq/${sample}

bowtie2 -x ${rRNA_index} \
-1 ${output_path}/quality_control/${sample}/${sample}.clean.1.fastq.gz \
-2 ${output_path}/quality_control/${sample}/${sample}.clean.2.fastq.gz \
--un-conc-gz ${output_path}/remove_rRNA/fastq/${sample}/${sample}.rm_rRNA.fq.gz \
-p 32 -S ${output_path}/remove_rRNA/fastq/${sample}/${sample}.aligned_rRNA.sam

rm ${output_path}/remove_rRNA/fastq/${sample}/${sample}.aligned_rRNA.sam
rm ${output_path}/quality_control/${sample}/${sample}.clean.1.fastq.gz
rm ${output_path}/quality_control/${sample}/${sample}.clean.2.fastq.gz


echo finish remove rRNA ${sample} `date`

# 4. mapping
echo start mapping ${sample} `date`
cd ${output_path}
mkdir -p mapping_expression/${sample}
cd mapping_expression/${sample}

STAR \
--runThreadN 32 \
--limitBAMsortRAM 20000000000 \
--outFilterType BySJout \
--outFilterMismatchNmax 10  \
--genomeDir ${genome_index_path} \
--readFilesIn ${output_path}/remove_rRNA/fastq/${sample}/${sample}.rm_rRNA.fq.1.gz \
${output_path}/remove_rRNA/fastq/${sample}/${sample}.rm_rRNA.fq.2.gz \
--readFilesCommand 'zcat' \
--outFileNamePrefix  ${sample} \
--outSAMtype BAM Unsorted \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMattributes All  --outSAMstrandField intronMotif --outBAMcompression 6 --outReadsUnmapped Fastx

samtools sort -T \
${output_path}/mapping_expression/${sample}/${sample}Aligned.out.sorted \
-o ${output_path}/mapping_expression/${sample}/${sample}Aligned.sortedByCoord.out.bam \
${output_path}/mapping_expression/${sample}/${sample}Aligned.out.bam

samtools sort -T \
${output_path}/mapping_expression/${sample}/${sample}Aligned.toTranscriptome.out.sorted \
-o ${output_path}/mapping_expression/${sample}/${sample}Aligned.toTranscriptome.out.sorted.bam \
${output_path}/mapping_expression/${sample}/${sample}Aligned.toTranscriptome.out.bam

samtools index \
${output_path}/mapping_expression/${sample}/${sample}Aligned.sortedByCoord.out.bam

samtools index \
${output_path}/mapping_expression/${sample}/${sample}Aligned.toTranscriptome.out.sorted.bam

rm ${output_path}/remove_rRNA/fastq/${sample}/${sample}.rm_rRNA.fq.1.gz
rm ${output_path}/remove_rRNA/fastq/${sample}/${sample}.rm_rRNA.fq.2.gz

echo finish mapping ${sample} `date`

# 5. featureCounts
echo start featureCounts ${sample}  `date`

cd ${output_path}
mkdir -p read_counts/${sample}
cd read_counts/${sample}

featureCounts \
-T 32 \
-s 0 \
-p -t CDS \
-g gene_id \
-a ${GTF} \
-o ${output_path}/read_counts/${sample}/${sample}.featurecounts.txt \
${output_path}/mapping_expression/${sample}/${sample}Aligned.sortedByCoord.out.bam

featureCounts \
-T 32 \
-s 0 \
-p -t exon \
-g gene_id \
-a ${GTF} \
-o ${output_path}/read_counts/${sample}/${sample}.featurecounts.all.txt \
${output_path}/mapping_expression/${sample}/${sample}Aligned.sortedByCoord.out.bam

echo finish featureCounts ${sample}  `date`

# 6.merge
cd ${output_path}
mkdir -p read_counts/result
cd read_counts/result

echo -e "gene_id	${sample}" >${output_path}/read_counts/result/${sample}.txt
cat ${output_path}/read_counts/${sample}/${sample}.featurecounts.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> ${output_path}/read_counts/result/${sample}.txt

echo -e "gene_id    ${sample}" >${output_path}/read_counts/result/${sample}.all.txt
cat ${output_path}/read_counts/${sample}/${sample}.featurecounts.all.txt| grep -v '#' | grep -v 'Geneid' | cut -f 1,7 >> ${output_path}/read_counts/result/${sample}.all.txt

echo finish ${sample} `date`

### 清除所有中间文件
#cd ${output_path}
#rm -r fastqc/${sample}
#rm -r quality_control/${sample}
#rm -r remove_rRNA/fastq/${sample}
#rm -r mapping_expression/${sample}
