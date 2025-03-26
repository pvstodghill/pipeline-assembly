import os
import glob

configfile: "config.yaml"

DATA=config['data'] if 'data' in config else "data"
BIN=DATA+"/bin"
PIPELINE=os.path.dirname(workflow.snakefile)

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

RAW_ASSEMBLY_FASTA = DATA+"/autocycler/consensus_assembly.fasta" if config['method'] == 'autocycler' else DATA+"/flye/assembly.fasta" if config['method'] == 'flye' else "error"
RAW_ASSEMBLY_GFA = DATA+"/autocycler/consensus_assembly.gfa" if config['method'] == 'autocycler' else DATA+"/flye/assembly_graph.gfa" if config['method'] == 'flye' else "error"

# ------------------------------------------------------------------------
# collect the inputs
# ------------------------------------------------------------------------

def flatten(alist):
    if alist == []:
        return []
    elif type(alist) is not list:
        return [alist]
    else:
        return flatten(alist[0]) + flatten(alist[1:])
    
def get_input_files(name):
    l = config[name]
    if not isinstance(l, list):
        l = [l]
    l = [os.path.expanduser(p) for p in l]
    l = [glob.glob(p) for p in l]
    l = flatten(l)
    if l == []:
        raise FileNotFoundError('input files for \''+name+'\' not found.')
    return l

rule make_raw_nanopore:
    input: get_input_files('nanopore')
    output: DATA+"/inputs/raw_nanopore.fastq.gz"
    shell: "cat {input} > {output}"

if get_config('short_R1') != None:
    rule make_raw_short_R1:
        input: get_input_files('short_R1')
        output: DATA+"/inputs/raw_short_R1.fastq.gz"
        shell: "cat {input} > {output}"

if get_config('short_R2') != None:
    rule make_raw_short_R2:
        input: get_input_files('short_R2')
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
# autocycler method
# ------------------------------------------------------------------------

# === download autocycler and auxiliary scripts ===

FLYE_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/flye.sh'
MINIASM_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/miniasm.sh'
RAVEN_SH_URL='https://raw.githubusercontent.com/rrwick/Autocycler/refs/heads/main/scripts/raven.sh'


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

# === autocycler subsample ===

rule make_subsamples:
    input:
        long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        gs=DATA+"/genome_size.txt"
    output: SUBSAMPLES_FQ
    params:
        subsamples = NUM_SUBSAMPLES
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler subsample \
        		  --reads {input.long_fq} \
			  --out_dir $(dirname {output[0]}) \
        		  --genome_size $(cat {input.gs}) \
        		  --count "{params.subsamples}"
        """

# === Run assembler {name} on subsample {i} ===

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

# === autocycler compress ===

rule run_autocycler_compress:
    input: 
        assemblies=ASSEMBLIES_FA
    output: DATA+"/autocycler/input_assemblies.gfa"
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler compress \
            --assemblies_dir $(dirname {input.assemblies[0]}) \
            --autocycler_dir $(dirname {output})
        """

# === autocycler cluster ===

checkpoint run_autocycler_cluster:
    input: 
        compresses_assemblies=DATA+"/autocycler/input_assemblies.gfa"
    output: directory(DATA+"/autocycler/clustering")
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler cluster \
            --autocycler_dir $(dirname {input.compresses_assemblies})
        """

def list_of_clusters(wildcards):
    ckp = checkpoints.run_autocycler_cluster.get(**wildcards).output[0]
    clusters = glob.glob(DATA+"/autocycler/clustering/qc_pass/cluster_*")
    return list(clusters)

# === autocycler trim ===

rule run_autocycler_trim:
    input:
        untrimmed="{cluster}/1_untrimmed.gfa"
    output: "{cluster}/2_trimmed.gfa"
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler trim --cluster_dir $(dirname {input.untrimmed})
        """

# === autocycler resolve ===

rule run_autocycler_resolve:
    input:
        trimmed="{cluster}/2_trimmed.gfa"
    output: "{cluster}/5_final.gfa"
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler resolve --cluster_dir $(dirname {input.trimmed})
        """

# === autocycler combine ===

rule run_autocycler_combine:
    input:
        gfas=expand("{cluster}/5_final.gfa",cluster=list_of_clusters)
    output:
        fasta=DATA+"/autocycler/consensus_assembly.fasta",
        gfa=DATA+"/autocycler/consensus_assembly.gfa"
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler combine \
            --autocycler_dir $(dirname {output.fasta}) \
            --in_gfas {input.gfas}
        """

# === generate version info ===

rule autocycler_version:
    output: DATA+"/versions/autocycler.txt"
    conda: "envs/autocycler.yaml"
    shell:
        """
        autocycler --version 2>&1 | tee {output}
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
# flye method
# ------------------------------------------------------------------------

rule run_flye:
    input:
        long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        gs=DATA+"/genome_size.txt"
    output: DATA+"/flye/assembly.fasta"
    params:
        args = get_config('flye_args','')
    threads: 9999
    conda: "envs/flye.yaml"
    shell:
        """
        flye \
            --nano-hq {input.long_fq} \
            --genome-size $(cat {input.gs}) \
            --out-dir $(dirname {output}) \
            --threads {threads} \
            {params.args}
        """

# ------------------------------------------------------------------------
# medaka
# ------------------------------------------------------------------------

rule run_medaka:
    input:
        reads=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        unpolished=RAW_ASSEMBLY_FASTA
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
            -d {input.unpolished} \
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

if (get_config('short_R1') != None) and (get_config('short_R2') != None):

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
        threads: 16 # fastp maxes at 16
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

if (get_config('short_R1') != None) and (get_config('short_R2') != None):

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

if (get_config('short_R1') != None) and (get_config('short_R2') != None):

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
# Create pipeline output assembly
# ------------------------------------------------------------------------

if (get_config('short_R1') != None) and (get_config('short_R2') != None):

    rule create_intermediate_fasta:
        input: DATA+"/pypolca/pypolca_corrected.fasta"
        output: DATA+"/intermediate.fasta"
        shell: "cp -a {input} {output}"

else:
    
    rule create_intermediate_fasta:
        input: DATA+"/medaka/consensus.fasta"
        output: DATA+"/intermediate.fasta"
        shell: "cp -a {input} {output}"

rule create_intermediate_gfa:
    input: {RAW_ASSEMBLY_GFA}
    output: DATA+"/intermediate.gfa"
    shell: "cp -a {input} {output}"

# ------------------------------------------------------------------------
# Run ReferenceSeeker
# ------------------------------------------------------------------------

if get_config('refseek_dir') != None:
    rule run_referenceseeker:
        input: DATA+"/intermediate.fasta"
        output: DATA+"/referenceseeker.log"
        params:
            refseek_dir=os.path.expanduser(get_config('refseek_dir')),
            refseek_args=get_config('refseek_args','-a 0.90 -u -r -g')
        conda: "envs/referenceseeker.yaml"
        threads: 9999
        shell:
            """
            REFSEEK={params.refseek_dir} \
            {PIPELINE}/scripts/run-referenceseeker -t {threads} {params.refseek_args} {input} \
                2>&1 | tee {output}
            """

rule referenceseeker_version:
    output: DATA+"/versions/referenceseeker.txt"
    conda: "envs/referenceseeker.yaml"
    shell:
        """
        referenceseeker --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# Run Unicycler
# ------------------------------------------------------------------------

rule run_unicycler:
    input:
        short_r1=DATA+"/fastp/trimmed_R1.fastq.gz",
        short_r2=DATA+"/fastp/trimmed_R2.fastq.gz",
        long_reads=DATA+"/filtlong/filtered_nanopore.fastq.gz",
    output: DATA+"/unicycler/assembly.fasta"
    threads: 9999
    conda: "envs/unicycler.yaml"
    shell:
        """
        unicycler -t {threads} \
          -1 {input.short_r1} \
          -2 {input.short_r2} \
          -l {input.long_reads} \
          -o $(dirname {output})
        """

rule unicycler_version:
    output: DATA+"/versions/unicycler.txt"
    conda: "envs/unicycler.yaml"
    shell:
        """
        unicycler --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# Compare input genome and Unicycler results with `dnadiff`
# ------------------------------------------------------------------------

rule run_dnadiff:
    input:
        raw=DATA+"/intermediate.fasta",
        unic=DATA+"/unicycler/assembly.fasta"
    output: DATA+"/dnadiff/out.report"
    conda: "envs/mummer4.yaml"
    shell:
        """
        dir=$(dirname {output})
        cp {input.raw} $dir/raw.fasta
        cp {input.unic} $dir/unicycler.fasta
        cd $dir
        dnadiff raw.fasta unicycler.fasta
        """

rule dnadiff_version:
    output: DATA+"/versions/dnadiff.txt"
    conda: "envs/mummer4.yaml"
    shell:
        """
        dnadiff --version 2>&1 | tee {output}
        """

# ------------------------------------------------------------------------
# Generate the summary
# ------------------------------------------------------------------------

rule make_summary:
    input:
        unpolished=RAW_ASSEMBLY_FASTA,
        intermediate_fasta=DATA+"/intermediate.fasta",
        intermediate_gfa=DATA+"/intermediate.gfa",
        referenceseeker_log=(DATA+"/referenceseeker.log" if 'refseek_dir' in config else []),
        autocycler_txt=DATA+"/versions/autocycler.txt" if config['method'] == 'autocycler' else [],
        fastp_txt=DATA+"/versions/fastp.txt" if 'short_R1' in config else [],
        filtlong_txt=DATA+"/versions/filtlong.txt",
        flye_txt=DATA+"/versions/flye.txt" if (config['method'] == 'autocycler' or config['method'] == 'flye') else [],
        lrge_txt=DATA+"/versions/lrge.txt" if 'genome_size' not in config else [],
        medaka_txt=DATA+"/versions/medaka.txt",
        miniasm_txt=DATA+"/versions/miniasm.txt" if config['method'] == 'autocycler' else [],
        polypolish_txt=DATA+"/versions/polypolish.txt" if 'short_R1' in config else [],
        pypolca_txt=DATA+"/versions/pypolca.txt" if 'short_R1' in config else [],
        raven_txt=DATA+"/versions/raven.txt" if config['method'] == 'autocycler' else [],
        unicycler_txt=DATA+"/versions/unicycler.txt",
        dnadiff_txt=DATA+"/versions/dnadiff.txt",
        referenceseeker_txt=DATA+"/versions/referenceseeker.txt" if 'refseek_dir' in config else [],
        dnadiff_report=(DATA+"/dnadiff/out.report" if 'skip_unicycler' not in config and 'short_R1' in config else []),
    output: DATA+"/summary-assembly.log"
    shell:
        """
        (
            echo
            echo === assembly summary ===
            if [ -e {DATA}/autocycler/consensus_assembly.fasta ] ; then
        	fgrep '>' {DATA}/autocycler/consensus_assembly.fasta
            elif [ -e {DATA}/flye/assembly_info.txt ] ; then
        	cat {DATA}/flye/assembly_info.txt
            else
        	fgrep '>' {input.unpolished}
            fi
            if [ "{input.referenceseeker_log}" ] ; then
                echo 
                echo === referenceseeker results ===
                fgrep -A10 '#ID' "{input.referenceseeker_log}"
            else
                echo 
                echo === no referenceseeker results ===
            fi
            if [ "{input.dnadiff_report}" ] ; then
                echo 
                echo === vs. unicycler ===
                head -n13 {input.dnadiff_report}  | tail -n+4
            fi
        ) | tee {output}
        """

# ------------------------------------------------------------------------
# Check Git status
# ------------------------------------------------------------------------

rule run_git:
    input: DATA+"/summary-assembly.log"
    output: DATA+"/git-assembly.log"
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
# Entry point
# ------------------------------------------------------------------------

rule all:
    input:
        DATA+"/summary-assembly.log",
        DATA+"/git-assembly.log"
    default_target: True

