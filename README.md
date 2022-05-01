# ont-assembly-snake

A snakemake-wrapper for creating *de novo* bacterial genome assemblies from Oxford Nanopore (ONT) sequencing data.
Combinations of assemblers and polishing tools can be applied to multiple samples at once.

Currently included programs:
* [flye](https://github.com/fenderglass/Flye)
* [raven](https://github.com/lbcb-sci/raven)
* [racon](https://github.com/lbcb-sci/racon)
* [medaka](https://github.com/nanoporetech/medaka)
* [pilon](https://github.com/broadinstitute/pilon/wiki)
* [Filtlong](https://github.com/rrwick/Filtlong)
* [homopolish](https://github.com/ythuang0522/homopolish)
* [proovframe](https://github.com/thackl/proovframe)

## Quick start
```bash
# Install
git clone https://github.com/pmenzel/ont-assembly-snake.git
conda config --add channels bioconda
conda env create -n ont-assembly-snake --file ont-assembly-snake/environment.yaml
conda activate ont-assembly-snake

# Prepare data
mkdir fastq-ont assemblies
cp /path/to/my/data/my_sample/ont_reads.fastq fastq-ont/mysample.fastq

# Declare desired combination of read filtering, assembly and polishing
mkdir assemblies/mysample_flye+racon2+medaka
mkdir assemblies/mysample+filtlong500_flye+racon2+medaka
mkdir assemblies/mysample_raven2+medaka

# Run workflow
snakemake -s ont-assembly-snake/Snakefile --cores 20
```


## Installation
Clone repository, for example into the existing folder `/opt/software/`:
```
git clone https://github.com/pmenzel/ont-assembly-snake.git /opt/software/ont-assembly-snake
```
Install [conda](https://docs.conda.io/en/latest/miniconda.html) and then create a new environment containing all the programs:
```
conda config --add channels bioconda
conda env create -n ont-assembly-snake --file /opt/software/ont-assembly-snake/environment.yaml
```
and activate environment:
```
conda activate ont-assembly-snake
```

## Usage
First, prepare a folder called `fastq-ont` containing the sequencing reads as
one fastq file per sample, e.g. `fastq-ont/sample1.fastq`.
For polishing with pilon using Illumina reads, the folder `fastq-illumina` must contain
those reads using `_R[12].fastq` suffixes, e.g. `fastq-illumina/sample1_R1.fastq` and `fastq-illumina/sample1_R2.fastq`.

Next, create a folder `assemblies` and inside create empty folders specifying
combinations of the desired combinations of read filtering, assembly and polishing steps.

Consecutive assembly and polishing steps need to be separated by a plus.

Read names and assembly+polishing need to be separated by an underscore.
NB: This also means that sample names must not contain underscores.

Both flye and raven do automatic polishing after assembly (by default 1 round
in flye and 2 rounds in raven). These defaults can be changed by
appending a number after the assembler name (e.g. `flye2`), where a 0 would
switch off the automatic polishing.

For running racon polishing multiple times, append a number specifying the
number of iterations, for example `sample1_flye+racon4` will run racon 4 times
on the flye assembly.

The ONT reads can be filtered by Filtlong prior to the assembly, for example
by reducing the input read set to contain only the highest quality reads up to a certain number of megabases.
This is specified by adding, i.e. for 500Mb, `+filtlong500` to the sample name.

The polishing steps with homopolish and proovframe always need to be the last one.

### Example folder structure
This example contains two samples with ONT sequencing reads and Illumina reads
for sample 2 only.  
We want to create various different assemblies for a later comparison, e.g.
with the [score-assemblies](https://github.com/pmenzel/score-assemblies)
snakemake pipeline.

For sample 1, the assembly should be done with flye (including the default single round of
polishing), followed by polishing the assembly with racon twice,
which in turn is again polished with medaka.
Sample 2 should be assembled by raven including two internal rounds polishing (raven uses racon internally),
followed by medaka and pilon (using the Illumina reads) polishing.
We also want to filter the ONT reads of sample 1 to only include the highest quality reads to a total of 500Mb
using Filtlong and apply the same assembly and polishing protocol.
```
.
├── fastq-ont
│   ├── sample1.fastq
│   └── sample2.fastq
│
├── fastq-illumina
│   ├── sample2_R1.fastq
│   └── sample2_R2.fastq
│
└── assemblies
    ├── sample1_flye+racon2+medaka
    ├── sample1+filtlong500_flye+racon2+medaka
    └── sample2_raven2+medaka+pilon
```
Snakemake will recursively handle the dependencies for each assembly, 
and create folders for all intermediate steps automatically.
For the above example, these folders will be created:
```
.
├── fastq-ont
│   ├── sample1.fastq
│   └── sample2.fastq
│
├── fastq-illumina
│   ├── sample2_R1.fastq
│   └── sample2_R2.fastq
│
└── assemblies
    ├── sample1_flye
    ├── sample1_flye+racon2
    ├── sample1_flye+racon2+medaka
    ├── sample1+filtlong500_flye
    ├── sample1+filtlong500_flye+racon2
    ├── sample1+filtlong500_flye+racon2+medaka
    ├── sample2_raven2
    ├── sample2_raven2+medaka
    └── sample2_raven2+medaka+pilon
```

### Run workflow

Run workflow in that folder, e.g. with 20 threads:
```
snakemake -k --use-conda -s /opt/software/ont-assembly-snake/Snakefile --cores 20
```

The assemblies are contained in the files `output.fa` in each folder.

Following additional configuration options can be passed to snakemake:

* `filtlong_min_read_length` (default: "1000"), used by Filtlong
* `medaka_model`, for specifying a medaka model to be used for all samples
* `map_medaka_model`, for specifying a tab-separated file with two columns mapping samples to medaka models

Options are specified using snakemake's `--config` parameter, for example:

```
snakemake -k --use-conda -s /opt/software/ont-assembly-snake/Snakefile --cores 20 --config map_medaka_model=map_medaka.tsv filtlong_min_read_length=5000
```
where `map_medaka.tsv` contains, for example, the two columns:
```
sample1     r941_min_high_g330
sample2     r941_min_high_g351
```

