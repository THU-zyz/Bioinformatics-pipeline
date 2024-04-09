shell.prefix('set -x; set -e;')

config = {"sample_id_path":"../Data/sample_ids.txt","indir":"../Data/","outdir":"output/test","sgRNA_library":"../Data/sgRNA_library_Metabolic.csv"}
sample_ids = open(config["sample_id_path"]).read().strip().split("\n")
indir = config["indir"]
outdir = config["outdir"]
sgRNA_library = config["sgRNA_library"]
# rule all 指定最终保留哪些文件
rule all:
    input:
        ""

rule fastqc:
    input:
        fastq1=indir+'/{sample_id}_1.fq.gz',
        fastq2=indir+'/{sample_id}_2.fq.gz'
    output:
        report1='{outdir}/fastqc/{sample_id}/{sample_id}_1_fastqc.html',
        report2='{outdir}/fastqc/{sample_id}/{sample_id}_2_fastqc.html'
    params:
        outdir='{outdir}/fastqc/{sample_id}'
    log:
        logfile='{outdir}/fastqc/{sample_id}/{sample_id}.log'
    shell:
        """
        fastqc {input.fastq1} \
        {input.fastq2} \
        --outdir {params.outdir} --noextract >> {log.logfile} 2>&1
        """

rule ungzip:
    input:
        fastq1=indir+'/{sample_id}_1.fq.gz',
        fastq2=indir+'/{sample_id}_2.fq.gz'
    output:
        fq1='{outdir}/ungzip/{sample_id}/{sample_id}_1.fq',
        fq2='{outdir}/ungzip/{sample_id}/{sample_id}_2.fq'

    shell:
        """
        gunzip -c {input.fastq1} > {output.fq1}
        gunzip -c {input.fastq2} > {output.fq2}
        """
rule pandaseq:
    input:
        fq1='{outdir}/ungzip/{sample_id}/{sample_id}_1.fq',
        fq2='{outdir}/ungzip/{sample_id}/{sample_id}_2.fq'
    output:
        merge_fq='{outdir}/merge/{sample_id}/{sample_id}.fq'
    shell:
    """
        pandaseq -F -f {input.fq1} \
        -r {input.fq2} \
        -w {output.merge_fq}
    """

rule mageck_count:
    input:
        sgRNA_library='{sgRNA_library}',
        merge_fq='{outdir}/merge/{sample_id}/{sample_id}.fq'
    output:
        mageck_count='{outdir}/mageck_count/{sample_id}/{sample_id}.count.txt'
    params:
        output_dir="{outdir}/mageck_count/{sample_id}/",
        name="{sample_id}",
        sample_label="{sample_id}"
    shell:
    """
        mageck count -l {input.sgRNA_library} \
        -n {params.name} --sample-label {params.sample_label} \
        --fastq {merge_fq} \
        
    """

    mageck count -l /public3/home/scg8403/Data/sgRNA/Data/sgRNA_library_Metabolic.csv -n PC9_all --sample-label "Treat1,Treat2,Treat3,Ctrl" --fastq /public3/home/scg8403/Data/sgRNA/Data/PG-1_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PG-2_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PG-3_L1.fq /public3/home/scg8403/Data/sgRNA/Data/PC_L1.fq --pdf-report
