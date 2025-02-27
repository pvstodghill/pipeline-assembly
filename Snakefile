import os
import glob

configfile: "config.yaml"

DATA=config['data'] if 'data' in config else "data"
BIN=DATA+"/bin"
PIPELINE=os.path.dirname(workflow.snakefile)

AUTOCYCLER_URL='https://github.com/rrwick/Autocycler/releases/download/v0.2.1/autocycler-linux-x86_64-gnu-v0.2.1.tar.gz'

FLYE_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/flye.sh'
MINIASM_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/miniasm.sh'
RAVEN_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/raven.sh'

GZIP="pigz"

def get_config(name, default=None):
    return config[name] if name in config else default

NUM_SUBSAMPLES=get_config("subsamples", "4")


SUBSAMPLES_IS=list(map((lambda j: "%02d" % (j)),range(1,int(NUM_SUBSAMPLES)+1)))

SUBSAMPLES_FQ = \
    expand(DATA+"/subsamples/sample_{i}.fastq",i=SUBSAMPLES_IS)

ASSEMBLERS=['flye','miniasm','raven']
ASSEMBLIES_FA = \
    expand(DATA+"/assemblies/{name}_{i}.fasta",
           name=ASSEMBLERS, i=SUBSAMPLES_IS)

# ------------------------------------------------------------------------
# download autocycler and auxiliary scripts
# ------------------------------------------------------------------------

rule download_autocycler:
    output: BIN+"/autocycler"
    shell:
        """
        cd $(dirname {output})
        wget -q {AUTOCYCLER_URL}
        tar xf autocycler*.tar.gz
        """

rule download_flye_script:
    output: BIN+"/flye.sh"
    shell:
        """
        wget -O {output} -q {FLYE_SH_URL}
        chmod +x {output}
        """

rule download_miniasm_script:
    output: BIN+"/miniasm.sh"
    shell:
        """
        wget -O {output} -q {MINIASM_SH_URL}
        chmod +x {output}
        """

rule download_raven_script:
    output: BIN+"/raven.sh"
    shell:
        """
        wget -O {output} -q {RAVEN_SH_URL}
        chmod +x {output}
        """

# ------------------------------------------------------------------------
# collect the inputs
# ------------------------------------------------------------------------

rule make_inputs_long_fq:
    input: os.path.expanduser(get_config("nanopore"))
    output: DATA+"/inputs/raw_nanopore.fastq.gz"
    shell:
        """
        cat {input} > {output}
        """

# ------------------------------------------------------------------------
# estimate the genome size
# ------------------------------------------------------------------------

if get_config('genome_size') != None:

    rule compute_genome_size:
        output: DATA+"/genome_size.txt"
        params: get_config('genome_size')
        shell:
            """
            echo {params} > {output}
            """

else:
    
    rule compute_genome_size:
        input: DATA+"/inputs/raw_nanopore.fastq.gz"
        output: DATA+"/genome_size.txt"
        conda: "envs/lrge.yaml"
        threads: 9999
        shell:
            """
            lrge -qqq --threads {threads} {input} > {output}
            """
            
# ------------------------------------------------------------------------
# run filtlong
# ------------------------------------------------------------------------


rule run_filtlong:
    input: DATA+"/inputs/raw_nanopore.fastq.gz",
    output: DATA+"/filtlong/filtered_nanopore.fastq.gz",
    threads: 9999
    conda: "envs/filtlong.yaml"
    shell:
        """
        filtlong --min_length 1000 --keep_percent 95 {input} \
            | {GZIP} > {output}
        """

# ------------------------------------------------------------------------
# autocycler subsample
# ------------------------------------------------------------------------

rule make_subsamples:
    input:
        autocycler=BIN+"/autocycler",
        long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        gs=DATA+"/genome_size.txt"
    output: SUBSAMPLES_FQ
    params:
        subsamples = NUM_SUBSAMPLES
    shell:
        """
        {input.autocycler} subsample \
        		  --reads {input.long_fq} \
			  --out_dir $(dirname {output[0]}) \
        		  --genome_size $(cat {input.gs}) \
        		  --count "{params.subsamples}"
        """

# ------------------------------------------------------------------------
# Run assembler {name} on subsample {i}
# ------------------------------------------------------------------------

rule make_one_assembly:
    input:
        script=BIN+"/{name}.sh",
        fq=DATA+"/subsamples/sample_{i}.fastq",
        gs=DATA+"/genome_size.txt"
    output: DATA+"/assemblies/{name}_{i}.fasta"
    threads: 9999
    conda: "envs/{name}.yaml"
    shell:
        """
        {input.script} \
            {input.fq} \
            $(dirname {output})/$(basename {output} .fasta) \
            {threads} $(cat {input.gs})
        """

# ------------------------------------------------------------------------
# autocycler compress
# ------------------------------------------------------------------------

rule run_autocycler_compress:
    input: 
        autocycler=BIN+"/autocycler",
        assemblies=ASSEMBLIES_FA
    output: DATA+"/autocycler/input_assemblies.gfa"
    shell:
        """
        {input.autocycler} compress \
            --assemblies_dir $(dirname {input.assemblies[0]}) \
            --autocycler_dir $(dirname {output})
        """

# ------------------------------------------------------------------------
# autocycler cluster
# ------------------------------------------------------------------------

checkpoint run_autocycler_cluster:
    input: 
        autocycler=BIN+"/autocycler",
        compresses_assemblies=DATA+"/autocycler/input_assemblies.gfa"
    output: directory(DATA+"/autocycler/clustering")
    shell:
        """
        {input.autocycler} cluster \
            --autocycler_dir $(dirname {input.compresses_assemblies})
        """

def list_of_clusters(wildcards):
    ckp = checkpoints.run_autocycler_cluster.get(**wildcards).output[0]
    clusters = glob.glob(DATA+"/autocycler/clustering/qc_pass/cluster_*")
    return list(clusters)

# ------------------------------------------------------------------------
# autocycler trim
# -----------------------------------------------------------------------

rule run_autocycler_trim:
    input:
        autocycler=BIN+"/autocycler",
        untrimmed="{cluster}/1_untrimmed.gfa"
    output: "{cluster}/2_trimmed.gfa"
    shell:
        """
        {input.autocycler} trim --cluster_dir $(dirname {input.untrimmed})
        """

# ------------------------------------------------------------------------
# autocycler resolve
# -----------------------------------------------------------------------

rule run_autocycler_resolve:
    input:
        autocycler=BIN+"/autocycler",
        trimmed="{cluster}/2_trimmed.gfa"
    output: "{cluster}/5_final.gfa"
    shell:
        """
        {input.autocycler} resolve --cluster_dir $(dirname {input.trimmed})
        """

# ------------------------------------------------------------------------
# autocycler combine
# ------------------------------------------------------------------------

rule run_autocycler_combine:
    input:
        autocycler=BIN+"/autocycler",
        gfas=expand("{cluster}/5_final.gfa",cluster=list_of_clusters)
    output: DATA+"/autocycler/consensus_assembly.fasta"
    shell:
        """
        {input.autocycler} combine \
            --autocycler_dir $(dirname {output}) \
            --in_gfas {input.gfas}
        """

# ------------------------------------------------------------------------
# medaka
# ------------------------------------------------------------------------

rule run_medaka:
    input:
        reads=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        consensus=DATA+"/autocycler/consensus_assembly.fasta"
    output: DATA+"/medaka/consensus.fasta"
    params:
        medaka_args="-b 10"
    threads: 16
    conda: "envs/medaka-gpu.yaml"
    shell:
        """
        medaka_consensus \
            {params.medaka_args} \
            -i {input.reads} \
            -d {input.consensus} \
            -o $(dirname {output}) \
            -t {threads} \
            --bacteria
        """

# ------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------

rule all:
    input: DATA+"/medaka/consensus.fasta"
    default_target: True

