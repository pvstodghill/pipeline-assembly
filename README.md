# pipeline-assembly

A pipeline for assembling bacteria genomes using one of several
methods.

The methods include,

- [Autocycler](https://github.com/rrwick/Autocycler) - fully
  automatic, consensus hybrid assembly.

## Cloning the repo

This pipeline uses Git submodules. The easiest way to clone this repo
(with a recent version of `git`) is

```
git clone --recurse-submodules https://github.com/pvstodghill/pipeline-autocycler.git
```

## Installing prereqs

You will want to install the following:

- [Snakemake](https://snakemake.github.io/)
- [Conda](https://conda.io)
- [Perl](https://www.perl.org/)

## Running on the example data

1. `( cd example ; bash download.bash )`

1. `snakemake --use-conda --configfile example/config.yaml`

## Running on your own data,

1. Copy `config.template.yaml` to `config.yaml`.  Edit `config.yaml`
   according to your needs and local environment.

1. `snakemake --use-conda`

