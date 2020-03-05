# ont-assembly-snake

A snakemake-wrapper for creating *de novo* bacterial genome assemblies from Oxford Nanopore (ONT)sequencing data.
Combinations of assemblers and polishing tools can be applied to multiple samples at once.

Currently included programs:
* [Flye](https://github.com/fenderglass/Flye)
* [raven](https://github.com/lbcb-sci/raven)
* [racon](https://github.com/lbcb-sci/racon)
* [medaka](https://github.com/nanoporetech/medaka)
* [pilon](https://github.com/broadinstitute/pilon/wiki)

## Installation
Clone repository, for example:
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
source activate ont-assembly-snake
```

## Usage
First, prepare a folder called `fastq-ont` containing the sequencing reads as
one fastq file per sample.
For polishing with Pilon using Illumina reads, the folder `fastq-illumina` must contain
those reads using `_R[12].fastq` suffixes.

Next, create a folder `assemblies` and inside create empty folders specifying
combinations of the desired assembly and polishing steps.

Consecutive steps need to be separated by a plus.

Sample names and assembly+polishing need to be separated by an underscore.
NB: This also means that sample names must not contain underscores.

For running racon polishing multiple times, append a number specifying the
number of iterations, for example `sample1_flye+racon4` will run racon 4 times
on the flye assembly.

### Example folder structure
This example contains two samples with ONT sequencing reads and Illumina reads
for sample2 only.  
We want to create various different assemblies for a later comparison, e.g.
with the [score-assemblies](https://github.com/pmenzel/score-assemblies)
snakemake pipeline.

For sample 1, the assembly should be done with flye, following by polishing the
assembly with racon once and 4 times, whereas the latter also being polished
with medaka. So we end up with 4 different assemblies for sample1.

Sample 2 should be assembled by raven including 2 internal steps of racon polishing,
followed by medaka as well as medaka and pilon (using the Illumina reads) polishing.
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
    ├── sample1_flye+racon
    ├── sample1_flye+racon4
    ├── sample1_flye+racon4+medaka
    ├── sample2_raven2
    ├── sample2_raven2+medaka
    └── sample2_raven2+medaka+pilon
```

Since Snakemake will handle the dependencies for each step, it is enough to
declare the final desired combination of assembly+polishing and all
folders for all intermediate steps will be created automatically.

Therefore, it is enough to manually create the following folders for the above example:
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
    ├── sample1_flye+racon
    ├── sample1_flye+racon4+medaka
    └── sample2_raven2+medaka+pilon
```


### Run workflow

Run workflow in that folder, e.g. with 20 threads:
```
snakemake -k -s /opt/software/ont-assembly-snake/Snakefile --cores 20
```

The assemblies are contained in the files `output.fa` in each folder.


Following additional configuration options can be passed to snakemake:

* `genome_size` (default: "6m"), used by flye
* `flye_iterations` (default: 4), number of polishing iterations in flye
* `medaka_model` (default: "r941_min_high_g344"), used by medaka

Options are specified using snakemake's `--config` parameter, for example:

```
snakemake -k -s /opt/software/ont-assembly-snake/Snakefile --cores 20 --config genome_size=5m
```

