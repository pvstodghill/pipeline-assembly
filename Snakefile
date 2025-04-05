import os
import glob

configfile: "config.yaml"

DATA=config['data'] if 'data' in config else "data"
BIN=DATA+"/bin"
PIPELINE=os.path.dirname(workflow.snakefile)

GZIP="pigz"

def get_config(name, default=None):
    return config[name] if name in config else default

NUM_SUBSAMPLES=int(get_config("subsamples", "4"))

ASSEMBLERS=['flye','miniasm','raven']

def compute_raw_assembly_fasta():
    if config['method'] == 'autocycler':
        return DATA+"/autocycler/consensus_assembly.fasta"
    elif config['method'] == 'trycycler':
        return DATA+"/trycycler/consensus.fasta"
    elif config['method'] == 'flye':
        return DATA+"/flye/assembly.fasta"
    else:
        return []

RAW_ASSEMBLY_FASTA = compute_raw_assembly_fasta()

def compute_raw_assembly_gfa():
    if config['method'] == 'autocycler':
        return DATA+"/autocycler/consensus_assembly.gfa"
    elif config['method'] == 'flye':
        return DATA+"/flye/assembly_graph.gfa"
    elif 'input_gfa' in config:
        return DATA+"/inputs/assembly.gfa"
    else:
        return []

#RAW_ASSEMBLY_GFA = compute_raw_assembly_gfa()

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

if config['method'] == 'external':
    rule make_raw_input_assembly_fasta:
        input: get_input_files('input_fasta')
        output: DATA+"/inputs/assembly.fasta"
        shell: "cat {input} > {output}"

if config['method'] == 'external' and 'input_gfa' in config:
    rule make_raw_input_assembly_gfa:
        input: get_input_files('input_gfa')
        output: DATA+"/inputs/assembly.gfa"
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
# method = flye
# ------------------------------------------------------------------------

rule run_flye:
    input:
        long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
        gs=DATA+"/genome_size.txt"
    output:
        fasta=DATA+"/flye/assembly.fasta",
        gs=DATA+"/flye/assembly_graph.gfa"
    params:
        args = get_config('flye_args','')
    threads: 9999
    conda: "envs/flye.yaml"
    shell:
        """
        flye \
            --nano-hq {input.long_fq} \
            --genome-size $(cat {input.gs}) \
            --out-dir $(dirname {output.fasta}) \
            --threads {threads} \
            {params.args}
        """

# ------------------------------------------------------------------------
# assembly scripts
# ------------------------------------------------------------------------

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

# ------------------------------------------------------------------------
# subsample the long reads
# ------------------------------------------------------------------------

if config['method'] == 'autocycler':

    rule make_subsamples:
        input:
            long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
            gs=DATA+"/genome_size.txt"
        output: expand(DATA+"/subsamples/sample_{i}.fastq",i=list(map((lambda j: "%02d" % (j)),range(1,NUM_SUBSAMPLES+1))))
        params:
            subsamples = NUM_SUBSAMPLES,
            min_read_depth = get_config('min_read_depth',25)
        conda: "envs/autocycler.yaml"
        shell:
            """
            autocycler subsample \
                              --reads {input.long_fq} \
                              --out_dir $(dirname {output[0]}) \
                              --genome_size $(cat {input.gs}) \
                              --count {params.subsamples} \
	       		      --min_read_depth {params.min_read_depth}
            """

elif config['method'] == 'trycycler':

    rule make_subsamples:
        input:
            long_fq=DATA+"/filtlong/filtered_nanopore.fastq.gz",
            gs=DATA+"/genome_size.txt"
        output: expand(DATA+"/subsamples/sample_{i}.fastq",i=list(map((lambda j: "%02d" % (j)),range(1,NUM_SUBSAMPLES*len(ASSEMBLERS)+1))))
        params:
            subsamples = NUM_SUBSAMPLES*len(ASSEMBLERS),
            min_read_depth = get_config('min_read_depth',25)
        threads: 9999
        conda: "envs/trycycler.yaml"
        shell:
            """
            trycycler subsample \
                      --reads {input.long_fq} \
            	      --out_dir $(dirname {output[0]}) \
                      --genome_size $(cat {input.gs}) \
                      --count {params.subsamples}  \
	       	      --min_read_depth {params.min_read_depth} \
                      --threads {threads}
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
# method = autocycler
# ------------------------------------------------------------------------

# === autocycler compress ===

AUTO_ASSEMBLIES_FA = \
    expand(DATA+"/assemblies/{name}_{i}.fasta",
           name=ASSEMBLERS, i=list(map((lambda j: "%02d" % (j)),range(1,NUM_SUBSAMPLES+1))))

rule run_autocycler_compress:
    input: 
        assemblies=AUTO_ASSEMBLIES_FA
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

def list_of_autocycler_clusters(wildcards):
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

#rule run_autocycler_combine:
rule run_autocycler:
    input:
        gfas=expand("{cluster}/5_final.gfa",cluster=list_of_autocycler_clusters)
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
# method = trycycler
# ------------------------------------------------------------------------

TRY_ASSEMBLIES_FA= \
        expand(DATA+"/assemblies/flye_{j}.fasta", j=list(map((lambda j: "%02d" % (j)),range(1,1*NUM_SUBSAMPLES+1)))) \
        + expand(DATA+"/assemblies/miniasm_{j}.fasta", j=list(map((lambda j: "%02d" % (j)),range(1*NUM_SUBSAMPLES+1,2*NUM_SUBSAMPLES+1)))) \
        + expand(DATA+"/assemblies/raven_{j}.fasta", j=list(map((lambda j: "%02d" % (j)),range(2*NUM_SUBSAMPLES+1,3*NUM_SUBSAMPLES+1))))

# === cluster the contigs from the different assembles ===
rule run_trycycler_cluster:
    input:
        assemblies=TRY_ASSEMBLIES_FA,
        long_reads=DATA+"/filtlong/filtered_nanopore.fastq.gz"
    output: directory(DATA+"/trycycler/raw_clusters")
    params: config['cluster_args'] if 'cluster_args' in config else ''
    threads: 9999
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler cluster \
            {params} \
            --threads {threads} \
            --assemblies {input.assemblies} \
            --reads {input.long_reads} \
            --out_dir {output}
        """
        
    
# === copy clusters; eliminate clusters, contigs, etc., as needed ===

# Cribbed from <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution>

checkpoint make_trycycler_clusters:
    input:
        clusters=DATA+"/trycycler/raw_clusters"
    output: DATA+"/trycycler/clusters/.done.txt"
    params:
        clusters_to_remove = get_config('remove_clusters',''),
        contigs_to_remove = get_config('remove_contigs',''),
        assemblies_to_remove = get_config('remove_assemblies',''),
    shell:
        """
        dir=$(dirname {output})
        rm -rf $dir
        cp --archive {input.clusters} $dir

        # remove clusters
        for index in {params.clusters_to_remove} ; do
            if [ -e $dir/${{index}} ] ; then
                path=$dir/${{index}}
            elif [ -e $dir/cluster_00${{index}} ] ; then
                path=$dir/cluster_00${{index}}
            elif [ -e $dir/cluster_0${{index}} ] ; then
                path=$dir/cluster_0${{index}}
            elif [ -e $dir/cluster_${{index}} ] ; then
                path=$dir/cluster_${{index}}
            fi
            if [ "$path" ] ; then
                echo "## removing cluster $index"
                (
                    set -x
                    rm -rf $path
                )
            fi
        done

        # remove contigs
        (
            shopt -s nullglob
            for name in {params.contigs_to_remove} ; do
                echo "## removing contig $name"
                paths="$(echo $dir/cluster_*/1_contigs/$name.fasta)"
                for path in $paths ; do
                    (
                        set -x
                        rm -f $path
                    )
                    cluster_xxx=$(dirname $(dirname $path))
                    rm -f ${{cluster_xxx}}/2_all_seqs.fasta
                done
            done
        )

        # remove assemblies
        (
            shopt -s nullglob
            for letter in {params.assemblies_to_remove} ; do
                echo "## removing assembly $letter"
                paths="$(echo $dir/cluster_*/1_contigs/${{letter}}_*.fasta)"
                for path in $paths ; do
                    (
                        set -x
                        rm -f $path
                    )
                    cluster_xxx=$(dirname $(dirname $path))
                    rm -f ${{cluster_xxx}}/2_all_seqs.fasta
                done
            done
        )

        touch {output}
        """

def list_of_trycycler_clusters(wildcards):
    ckp = checkpoints.make_trycycler_clusters.get(**wildcards).output[0]
    clusters = glob.glob(DATA+"/trycycler/clusters/cluster_*")
    return list(clusters)
    #return expand(ckp_reconciled_dir+"/cluster_{i}",i=clusters)

# === reconcile the contigs in each cluster ===

rule run_trycycler_reconcile:
    input:
        contigs=DATA+"/trycycler/clusters/{cluster}/1_contigs",
        long_reads=DATA+"/filtlong/filtered_nanopore.fastq.gz"
    output: DATA+"/trycycler/clusters/{cluster}/2_all_seqs.fasta"
    params:
        args_global=get_config('reconcile_args',''),
        args_cluster=get_config('reconcile_args_{cluster}','')
    threads: 9999
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler reconcile \
                --threads {threads} \
                --reads {input.long_reads} \
                --cluster_dir $(dirname {input.contigs}) \
                {params.args_global} {params.args_cluster}
        """

# === multiple sequence alignment (MSA) ===

rule run_trycycler_msa:
    input: DATA+"/trycycler/clusters/{cluster}/2_all_seqs.fasta"
    output: DATA+"/trycycler/clusters/{cluster}/3_msa.fasta"
    threads: 9999
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler msa \
                --threads {threads} \
                --cluster_dir $(dirname {input})
        """

# === Partition long reads across clusters ===

rule run_trycycler_partition:
    input:
        contigs=expand("{cluster}/3_msa.fasta",cluster=list_of_trycycler_clusters),
        long_reads=DATA+"/filtlong/filtered_nanopore.fastq.gz"
    output: DATA+"/trycycler/clusters/.done2.txt"
    threads: 9999
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler partition \
                 --threads {threads} \
                 --reads {input.long_reads} \
                 --cluster_dirs $(for c in {input.contigs} ; do dirname $c ; done)
        touch {output}
        """

# === Generate a consensus sequence for each cluster ===

rule run_trycycler_consensus:
    input:
        DATA+"/trycycler/clusters/.done2.txt",
        #DATA+"/trycycler/clusters/{cluster}/4_reads.fastq"
    output:
        gfa=DATA+"/trycycler/clusters/{cluster}/5_chunked_sequence.gfa",
        tmp_fna=DATA+"/trycycler/clusters/{cluster}/6_initial_consensus.fasta",
        fna=DATA+"/trycycler/clusters/{cluster}/7_final_consensus.fasta",
    threads: 9999
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler consensus \
                 --threads {threads} \
                 --cluster_dir $(dirname {output.fna})
        """

# === Generate the aggregate consensus assembly ===

rule make_trycycler_consensus_fasta:
    input: expand("{cluster}/7_final_consensus.fasta",cluster=list_of_trycycler_clusters)
    output: DATA+"/trycycler/consensus.fasta"
    shell: "cat {input} > {output}"


# === generate version info ===

rule trycycler_version:
    output: DATA+"/versions/trycycler.txt"
    conda: "envs/trycycler.yaml"
    shell:
        """
        trycycler --version 2>&1 | tee {output}
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

if config['method'] == 'external':
    
    rule create_intermediate_fasta:
        input: DATA+"/inputs/assembly.fasta"
        output: DATA+"/intermediate.fasta"
        shell: "cp -a {input} {output}"

elif (get_config('short_R1') != None) and (get_config('short_R2') != None):

    rule create_intermediate_fasta:
        input: DATA+"/pypolca/pypolca_corrected.fasta"
        output: DATA+"/intermediate.fasta"
        shell: "cp -a {input} {output}"

else:
    
    rule create_intermediate_fasta:
        input: DATA+"/medaka/consensus.fasta"
        output: DATA+"/intermediate.fasta"
        shell: "cp -a {input} {output}"


if compute_raw_assembly_gfa() != []:

    rule create_intermediate_gfa:
        input: compute_raw_assembly_gfa()
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
        #unpolished=RAW_ASSEMBLY_FASTA,
        intermediate_fasta=DATA+"/intermediate.fasta",
        intermediate_gfa=DATA+"/intermediate.gfa" if compute_raw_assembly_gfa() != [] else [],
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
    params:
        method=config['method'],
    shell:
        """
        (
            echo
            echo === assembly summary ===
            case {params.method} in
            autocycler)
        	fgrep '>' {DATA}/autocycler/consensus_assembly.fasta
        	;;
            flye)
        	cat {DATA}/flye/assembly_info.txt
		;;
            external)
        	if [ "{input.intermediate_gfa}" ] ; then
        	    echo 1>&2 implementme: topology from gfa
        	    exit 1
        	else
        	    cat {input.intermediate_fasta} | {PIPELINE}/scripts/fasta-length -f
        	fi
        	;;
            trycycler)
        	cat {input.intermediate_fasta} | {PIPELINE}/scripts/fasta-length -f
        	;;
            *)
        	echo 1>&2 cannot happen: {params.method}
        	exit 1
            esac
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

