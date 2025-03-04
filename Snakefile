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

rule make_raw_nanopore:
    input: os.path.expanduser(config['nanopore'])
    output: DATA+"/inputs/raw_nanopore.fastq.gz"
    shell: "cat {input} > {output}"

rule make_raw_short_R1:
    input: os.path.expanduser(config['short_R1'])
    output: DATA+"/inputs/raw_short_R1.fastq.gz"
    shell: "cat {input} > {output}"

rule make_raw_short_R2:
    input: os.path.expanduser(config['short_R2'])
    output: DATA+"/inputs/raw_short_R2.fastq.gz"
    shell: "cat {input} > {output}"

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
    rule lrge_version:
        output: DATA+"/versions/lrge.txt"
        conda: "envs/lrge.yaml"
        shell:
            """
            lrge --version 2>&1 | tee {output}
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

rule filtlong_version:
    output: DATA+"/versions/filtlong.txt"
    conda: "envs/filtlong.yaml"
    shell:
        """
        filtlong --version 2>&1 | tee {output}
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

rule autocycler_version:
    input:
        autocycler=BIN+"/autocycler",
    output: DATA+"/versions/autocycler.txt"
    shell:
        """
        {input.autocycler} --version 2>&1 | tee {output}
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

rule flye_version:
    output: DATA+"/versions/flye.txt"
    conda: "envs/flye.yaml"
    shell:
        """
        flye --version 2>&1 | tee {output}
        """

rule miniasm_version:
    output: DATA+"/versions/miniasm.txt"
    conda: "envs/miniasm.yaml"
    shell:
        """
        miniasm -V 2>&1 | tee {output}
        """

rule raven_version:
    output: DATA+"/versions/raven.txt"
    conda: "envs/raven.yaml"
    shell:
        """
        raven --version 2>&1 | tee {output}
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

rule medaka_version:
    output: DATA+"/versions/medaka.txt"
    conda: "envs/medaka-gpu.yaml"
    shell:
        """
        medaka --version 2>&1 | tee {output}
        """
# ------------------------------------------------------------------------
# fastp
# ------------------------------------------------------------------------

rule run_fastp:
    input:
        r1=DATA+"/inputs/raw_short_R1.fastq.gz",
        r2=DATA+"/inputs/raw_short_R2.fastq.gz"
    output:
        r1=DATA+"/fastp/trimmed_R1.fastq.gz",
        r2=DATA+"/fastp/trimmed_R2.fastq.gz",
        u=DATA+"/fastp/u.fastq.gz",
        json=DATA+"/fastp/fastp.json",
        html=DATA+"/fastp/fastp.html",
    threads: 9999
    conda: "envs/fastp.yaml"
    shell:
        """
        fastp \
            --thread {threads} \
            --json {output.json} --html {output.html} \
            --in1 {input.r1} --in2 {input.r2}  \
            --out1 {output.r1} --out2 {output.r2} \
            --unpaired1 {output.u} --unpaired2 {output.u}
        """

rule fastp_version:
    output: DATA+"/versions/fastp.txt"
    conda: "envs/fastp.yaml"
    shell:
        """
        fastp --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# polypolish
# ------------------------------------------------------------------------

rule run_polypolish:
    input:
        draft=DATA+"/medaka/consensus.fasta",
        r1=DATA+"/fastp/trimmed_R1.fastq.gz",
        r2=DATA+"/fastp/trimmed_R2.fastq.gz",
    output: DATA+"/polypolish/polished.fasta"
    threads: 9999
    conda: "envs/polypolish.yaml"
    shell:
        """
        dir=$(dirname {output})
        cp {input.draft} $dir/draft.fasta
        cp {input.r1} $dir/r1.fastq.gz
        cp {input.r2} $dir/r2.fastq.gz
        coverage=$({PIPELINE}/scripts/fastq-coverage -p 0 -g $dir/draft.fasta $dir/r1.fastq.gz $dir/r2.fastq.gz)
        cd $dir
        bwa index draft.fasta
        bwa mem -t {threads} -a draft.fasta r1.fastq.gz > align1.sam
        bwa mem -t {threads} -a draft.fasta r2.fastq.gz > align2.sam
        polypolish filter \
            --in1 align1.sam --in2 align2.sam \
            --out1 filtered1.sam --out2 filtered2.sam
        if [ $coverage -le 25 ] ; then
            CAREFUL=--careful
        else
            CAREFUL=
        fi
        polypolish polish $CAREFUL draft.fasta filtered1.sam filtered2.sam > polished.fasta
        """

rule polypolish_version:
    output: DATA+"/versions/polypolish.txt"
    conda: "envs/polypolish.yaml"
    shell:
        """
        polypolish --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# pypolca
# ------------------------------------------------------------------------

rule run_pypolca:
    input:
        draft=DATA+"/polypolish/polished.fasta",
        r1=DATA+"/fastp/trimmed_R1.fastq.gz",
        r2=DATA+"/fastp/trimmed_R2.fastq.gz",
    output: DATA+"/pypolca/pypolca_corrected.fasta"
    threads: 16 # weird samtools memory problem
    conda: "envs/pypolca.yaml"
    shell:
        """
        pypolca run --force -a {input.draft} -1 {input.r1} -2 {input.r2} \
            -t {threads} -o $(dirname {output}) --careful
        """

rule pypolca_version:
    output: DATA+"/versions/pypolca.txt"
    conda: "envs/pypolca.yaml"
    shell:
        """
        pypolca --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# Run ReferenceSeeker
# ------------------------------------------------------------------------

REFSEEK_CUTOFF = get_config('refseek_cutoff','0.95')

if get_config('refseek_dir') != None:
    rule run_referenceseeker:
        input: DATA+"/pypolca/pypolca_corrected.fasta"
        output: DATA+"/referenceseeker.log"
        params:
            refseek_dir=os.path.expanduser(get_config('refseek_dir')),
            refseek_cutoff="-u" if REFSEEK_CUTOFF == 0.0 else "-a "+REFSEEK_CUTOFF,
            refseek_dbs=get_config('refseek_dbs','-r')
        conda: "envs/referenceseeker.yaml"
        shell:
            """
            REFSEEK={params.refseek_dir} \
            {PIPELINE}/scripts/run-referenceseeker {params.refseek_cutoff} \
            		{params.refseq_dbs} {input} \
                | tee {output}
            """

rule referenceseeker_version:
    output: DATA+"/versions/referenceseeker.txt"
    conda: "envs/referenceseeker.yaml"
    shell:
        """
        referenceseeker --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# Check Git status
# ------------------------------------------------------------------------

rule run_git:
    input:
        DATA+"/versions/autocycler.txt",
        DATA+"/versions/fastp.txt",
        DATA+"/versions/filtlong.txt",
        DATA+"/versions/flye.txt",
        DATA+"/versions/lrge.txt",
        DATA+"/versions/medaka.txt",
        DATA+"/versions/miniasm.txt",
        DATA+"/versions/polypolish.txt",
        DATA+"/versions/pypolca.txt",
        DATA+"/versions/raven.txt",
        DATA+"/versions/referenceseeker.txt",
        DATA+"/pypolca/pypolca_corrected.fasta",
        (DATA+"/referenceseeker.log" if 'refseek_dir' in config else [])
    output: DATA+"/git-autocycler.log"
    shell:
        """
	(
	    cd {PIPELINE}
	    echo
	    ( set -x ; git status )
	    echo
	    ( set -x ; git log -n1 )
	) 2>&1 | tee {output}
        """ 

# ------------------------------------------------------------------------
# Generate the summary
# ------------------------------------------------------------------------

rule make_summary:
    input:
        autocycler_assembly=DATA+"/autocycler/consensus_assembly.fasta",
        dnadiff_report=DATA+"/dnadiff/out.report",
        referenceseeker_log=DATA+"/referenceseeker.log"
    output: DATA+"/summary-autocycler.log"
    shell:
        """
        (
            echo
            echo === autocycler summary ===
            fgrep '>' {input.autocycler_assembly}
            echo 
            echo === autocycler vs. unicycler ===
            head -n8 {input.dnadiff_report}  | tail -n+4
            echo 
            echo === referenceseeker results ===
            cat {input.referenceseeker_log}
        ) | tee {output}
        """
                

# ------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------

rule all:
    input: DATA+"/summary-autocycler.log"
    default_target: True

