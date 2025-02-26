import os
import glob

configfile: "config.yaml"

DATA=config['data'] if 'data' in config else "data"
PIPELINE=os.path.dirname(workflow.snakefile)

AUTOCYCLER_URL='https://github.com/rrwick/Autocycler/releases/download/v0.2.1/autocycler-linux-x86_64-gnu-v0.2.1.tar.gz'

# ------------------------------------------------------------------------

rule all:
    input:
        DATA+"/autocycler"

rule download_autocycler:
    output: DATA+"/autocycler"
    shell:
        """
        cd $(dirname {output})
        wget -q {AUTOCYCLER_URL}
        tar xf autocycler*.tar.gz
        """
