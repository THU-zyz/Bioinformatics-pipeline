# `-x` for debugging and `-e` for error handing
# These options enables debugging mode and the script to exit immediately if any command it runs exits with a non-zero status.
shell.prefix('set -x; set -e;')

config = {"sample_id_path":"","in_dir":"","out_dir":""}

# read all sample ids from sample id path
sample_ids = open(config['sample_id_path']).read().strip().split("\n")
indir = config['indir']
outdir = config['outdir']
genome = config['']


# rule all is a special rule. it is used to define the end goal or output file of the workflow.
rule all:
    # define input files in `all` pipeline.
    input:
        # Rules for replacing files in a snakemake workflow, similar to *.txt in shell.
        expand("{outdir}//{sample_id}/{sample_id}_1_fastqc.html",outdir=outrdir,sample_id=sample_ids)
        