#! /bin/bash -e

rm -f autocycler-demo-dataset.tar
rm -f reads.fastq.gz
rm -f truth.fasta

set -x

wget -q https://github.com/rrwick/Autocycler/releases/download/v0.1.0/autocycler-demo-dataset.tar
tar xvf autocycler-demo-dataset.tar 
