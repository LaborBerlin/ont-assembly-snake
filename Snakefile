from glob import glob
import pandas as pd

from snakemake.utils import min_version

min_version("6.0")

if config.get("run_score_assemblies", False):

    module score_assemblies:
        snakefile:
            github("pmenzel/score-assemblies", path="Snakefile", branch="master")
            #  "score-assemblies/Snakefile"
        config:
            config

    use rule * from score_assemblies as score_assemblies_*

    ruleorder: proovframe_diamond_index > score_assemblies_diamond_ref_makedb


filtlong_min_read_length = "1000"
if config.get("filtlong_min_read_length", False):
    filtlong_min_read_length = config["filtlong_min_read_length"]
print("filtlong min. read length = " + str(filtlong_min_read_length))

medaka_model = None
if config.get("medaka_model", False):
    medaka_model = config["medaka_model"]
    print("Medaka model = " + medaka_model)

map_medaka_model = None
if config.get("map_medaka_model", False):
    medaka_model_file = config["map_medaka_model"]
    table = pd.read_csv(medaka_model_file, sep="\t", lineterminator="\n", header=None)
    map_medaka_model = dict(zip(table[0], table[1]))

# this option is needed by Canu
target_genome_size = None
if config.get("genome_size", False):
    target_genome_size = config["genome_size"]
    print("Target genome size = " + str(target_genome_size) + "M")


wildcard_constraints:
    sample="[^_]+",
    assembly="[^_]+",
    sample_assembly="[^/]+",
    num="[0-9]+",


def get_ont_fq(wildcards):
    if "filtlong" in wildcards.sample:
        return "fastq-ont/" + wildcards.sample + ".fastq"
    elif "rasusa" in wildcards.sample:
        return "fastq-ont/" + wildcards.sample + ".fastq"
    else:
        return glob("fastq-ont/" + wildcards.sample + ".fastq*")


# use split("+")[0] here for removing the +filtlong... or +rasusa... suffices from sample names for Illumina reads
def get_R1_fq(wildcards):
    return glob("fastq-illumina/" + wildcards.sample.split("+")[0] + "_R1.fastq*")


def get_R2_fq(wildcards):
    return glob("fastq-illumina/" + wildcards.sample.split("+")[0] + "_R2.fastq*")


(references,) = glob_wildcards("references/{ref,[^/\\\\]+}.fa")

(references_protein,) = glob_wildcards("references-protein/{ref,[^/\\\\]+}.faa")

(sample_assemblies,) = glob_wildcards("assemblies/{sample_assembly,[^/]+}/")
# ignore symlinks in assemblies/folder, e.g. sample_flye.fa -> assemblies/sample_flye/output.fa
sample_assemblies = [a for a in sample_assemblies if not re.search("\.fa", a)]

# if config files with list of assemblies is given, then use this instead of folders in assemblies/
if config.get("assemblies", False):
    sample_assemblies = list(set(config["assemblies"]))

# if any desired assembly requires homopolish then at least one reference genome should be provided
if not references and [
    string for string in sample_assemblies if "homopolish" in string
]:
    quit(
        "Error: must provide at least one reference genome sequence when using homopolish"
    )

# if any desired assembly requires proovframe then at least one reference proteome should be provided
if not references_protein and [
    string for string in sample_assemblies if "proovframe" in string
]:
    quit(
        "Error: must provide at least one reference protein file when using proovframe"
    )

# if any desired assembly uses Canu then need to provide target genome size
if not target_genome_size and [
    string for string in sample_assemblies if "_canu" in string
]:
    quit(
        "Error: must provide target genome size when using Canu assembler, use option (e.g. for 5.2Mb): --config genome_size=5.2"
    )

list_outputs = expand(
    "assemblies/{sample_assembly}/output.fa", sample_assembly=sample_assemblies
)
list_outputs_links = expand(
    "assemblies/{sample_assembly}.fa", sample_assembly=sample_assemblies
)
# remove homopolish and proovframe assemblies from default list
list_outputs = [
    i for i in list_outputs if not re.search("homopolish|proovframe", i, re.IGNORECASE)
]
list_outputs_links = [
    i
    for i in list_outputs_links
    if not re.search("homopolish|proovframe", i, re.IGNORECASE)
]

# make lists for homopolish, one entry for each reference genome
list_outputs_homopolish = expand(
    "assemblies/{sample_assembly}/output_{ref}.fa",
    ref=references,
    sample_assembly=[
        i for i in sample_assemblies if re.search("homopolish$", i, re.IGNORECASE)
    ],
)
list_outputs_links_homopolish = expand(
    "assemblies/{sample_assembly}{ref}.fa",
    ref=references,
    sample_assembly=[
        i for i in sample_assemblies if re.search("homopolish$", i, re.IGNORECASE)
    ],
)

# make lists for proovframe, one entry for each reference protein file
list_outputs_proovframe = expand(
    "assemblies/{sample_assembly}/output_{ref}.fa",
    ref=references_protein,
    sample_assembly=[
        i for i in sample_assemblies if re.search("proovframe$", i, re.IGNORECASE)
    ],
)
list_outputs_links_proovframe = expand(
    "assemblies/{sample_assembly}{ref}.fa",
    ref=references_protein,
    sample_assembly=[
        i for i in sample_assemblies if re.search("proovframe$", i, re.IGNORECASE)
    ],
)


rule all:
    input:
        list_outputs,
        list_outputs_links,
        list_outputs_homopolish,
        list_outputs_links_homopolish,
        list_outputs_proovframe,
        list_outputs_links_proovframe,


rule filtlong:
    threads: 1
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+filtlong.fastq",
    log:
        "fastq-ont/{sample}_filtlong_log.txt",
    shell:
        """
        filtlong --min_length {filtlong_min_read_length} {input} > {output} 2>{log}
        """


rule rasusaMB:
    conda:
        "env/conda-rasusa.yaml"
    threads: 1
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+rasusaMB{num}.fastq",
    log:
        "fastq-ont/{sample}_rasusaMB{num}_log.txt",
    shell:
        """
         rasusa --bases {wildcards.num}m -i {input} -o {output} 2>{log}
        """


rule filtlongMB:
    threads: 1
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+filtlongMB{num}.fastq",
    log:
        "fastq-ont/{sample}_filtlongMB{num}_log.txt",
    shell:
        """
        filtlong --min_length {filtlong_min_read_length} -t {wildcards.num}000000 {input} > {output} 2>{log}
        """


# for keeping num PerCent of the bases
rule filtlongPC:
    threads: 1
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+filtlongPC{num}.fastq",
    log:
        "fastq-ont/{sample}_filtlongPC{num}_log.txt",
    shell:
        """
        filtlong --min_length {filtlong_min_read_length} --keep_percent {wildcards.num} {input} > {output} 2>{log}
        """


rule filtlongMBql:
    threads: 1
    wildcard_constraints:
        mb="[0-9]+",
        qweight="[0-9]+",
        lweight="[0-9]+",
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+filtlongMB{mb},{qweight},{lweight}.fastq",
    log:
        "fastq-ont/{sample}_filtlongMB{mb},{qweight},{lweight}_log.txt",
    shell:
        """
        filtlong --min_length {filtlong_min_read_length} --mean_q_weight {wildcards.qweight} --length_weight {wildcards.lweight}  -t {wildcards.mb}000000 {input} > {output} 2>{log}
        """


rule filtlongMBqln:
    threads: 1
    wildcard_constraints:
        mb="[0-9]+",
        readlen="[0-9]+",
        qweight="[0-9]+",
        lweight="[0-9]+",
    input:
        fq=get_ont_fq,
    output:
        "fastq-ont/{sample}+filtlongMB{mb},{qweight},{lweight},{readlen}.fastq",
    log:
        "fastq-ont/{sample}_filtlongMB{mb},{qweight},{lweight},{readlen}_log.txt",
    shell:
        """
        filtlong --min_length {wildcards.readlen} --mean_q_weight {wildcards.qweight} --length_weight {wildcards.lweight}  -t {wildcards.mb}000000 {input} > {output} 2>{log}
        """


rule miniasm:
    conda:
        "env/conda-miniasm.yaml"
    threads: 10
    input:
        fqont=get_ont_fq,
    output:
        fa="assemblies/{sample}_miniasm/output.fa",
        link="assemblies/{sample}_miniasm.fa",
    log:
        "assemblies/{sample}_miniasm/log.txt",
    shell:
        """
        minimap2 -x ava-ont -t {threads} {input} {input} > assemblies/{wildcards.sample}_miniasm/minimap2_overlap.paf 2>{log}
        miniasm -f {input} assemblies/{wildcards.sample}_miniasm/minimap2_overlap.paf > assemblies/{wildcards.sample}_miniasm/minimap2_miniasm.gfa 2>>{log}
        perl -lsane 'print ">$F[1]\n$F[2]" if $F[0] =~ /S/;' assemblies/{wildcards.sample}_miniasm/minimap2_miniasm.gfa > {output.fa} 2>>{log}
        ln -sr {output.fa} {output.link}
        """


rule unicycler:
    conda:
        "env/conda-unicycler.yaml"
    threads: 10
    input:
        fqont=get_ont_fq,
        fq1=get_R1_fq,
        fq2=get_R2_fq,
    output:
        fa="assemblies/{sample}_unicycler/output.fa",
        link="assemblies/{sample}_unicycler.fa",
    log:
        "assemblies/{sample}_unicycler/log.txt",
    shell:
        """
        # del spades folder if already exists (e.g. when workflow was canceled), so that it does not warn about it upon restart
        [ -d "assemblies/{wildcards.sample}_unicycler/spades_assembly" ] && rm -r "assemblies/{wildcards.sample}_unicycler/spades_assembly" >{log}
        unicycler -1 {input.fq1} -2 {input.fq2} -l {input.fqont} -t {threads} --keep 0 -o assemblies/{wildcards.sample}_unicycler/ >>{log} 2>&1
        cp assemblies/{wildcards.sample}_unicycler/assembly.fasta {output.fa} 2>>{log}
        ln -sr {output.fa} {output.link}
        """


# flye with default number of polishing rounds (=1 in flye v2.9)
rule flye:
    conda:
        "env/conda-flye.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_flye/output.fa",
        link="assemblies/{sample}_flye.fa",
    log:
        "assemblies/{sample}_flye/log.txt",
    shell:
        """
        flye --nano-raw {input.fq} -o assemblies/{wildcards.sample}_flye/ -t {threads} 2>{log}
        mv assemblies/{wildcards.sample}_flye/assembly.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


rule flyeX:
    conda:
        "env/conda-flye.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_flye{num}/output.fa",
        link="assemblies/{sample}_flye{num}.fa",
    log:
        "assemblies/{sample}_flye{num}/log.txt",
    shell:
        """
        flye --nano-raw {input.fq} -o assemblies/{wildcards.sample}_flye{wildcards.num}/ -t {threads} -i {wildcards.num} 2>{log}
        mv assemblies/{wildcards.sample}_flye{wildcards.num}/assembly.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


# flye for high quality ONT reads (from flye 2.9), with default number of polishing rounds (=1 in flye v2.9)
rule flyehq:
    conda:
        "env/conda-flye.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_flyehq/output.fa",
        link="assemblies/{sample}_flyehq.fa",
    log:
        "assemblies/{sample}_flyehq/log.txt",
    shell:
        """
        flye --nano-hq {input.fq} -o assemblies/{wildcards.sample}_flyehq/ -t {threads} 2>{log}
        mv assemblies/{wildcards.sample}_flyehq/assembly.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


rule flyehqX:
    conda:
        "env/conda-flye.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_flyehq{num}/output.fa",
        link="assemblies/{sample}_flyehq{num}.fa",
    log:
        "assemblies/{sample}_flyehq{num}/log.txt",
    shell:
        """
        flye --nano-hq {input.fq} -o assemblies/{wildcards.sample}_flyehq{wildcards.num}/ -t {threads} -i {wildcards.num} 2>{log}
        mv assemblies/{wildcards.sample}_flyehq{wildcards.num}/assembly.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


# for running raven with default number of racon-polishing rounds (=2 in raven v0.0.8)
rule raven:
    conda:
        "env/conda-raven.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_raven/output.fa",
        link="assemblies/{sample}_raven.fa",
    log:
        "assemblies/{sample}_raven/log.txt",
    shell:
        """
        raven --disable-checkpoints -t {threads} {input.fq} >{output.fa} 2>{log}
        ln -sr {output.fa} {output.link}
        """


# for running raven with racon polishing X times
rule ravenX:
    conda:
        "env/conda-raven.yaml"
    threads: 5
    input:
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_raven{num}/output.fa",
        link="assemblies/{sample}_raven{num}.fa",
    log:
        "assemblies/{sample}_raven{num}/log.txt",
    shell:
        """
        raven --disable-checkpoints -p {wildcards.num} -t {threads} {input.fq} >{output.fa} 2>{log}
        ln -sr {output.fa} {output.link}
        """


rule canu:
    conda:
        "env/conda-canu.yaml"
    threads: 10
    input:
        fqont=get_ont_fq,
    output:
        fa="assemblies/{sample}_canu/output.fa",
        link="assemblies/{sample}_canu.fa",
    log:
        "assemblies/{sample}_canu/log.txt",
    shell:
        """
        canu -nanopore -d assemblies/{wildcards.sample}_canu/ -p output useGrid=false maxThreads={threads} genomeSize={target_genome_size}m {input.fqont} >>{log} 2>&1
        cp assemblies/{wildcards.sample}_canu/output.contigs.fasta {output.fa} 2>>{log}
        ln -sr {output.fa} {output.link}
        """


# for running racon once
rule racon:
    conda:
        "env/conda-racon.yaml"
    threads: 5
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_{assembly}+racon/output.fa",
        link="assemblies/{sample}_{assembly}+racon.fa",
        sam=temp("assemblies/{sample}_{assembly}+racon/map.sam"),
    log:
        "assemblies/{sample}_{assembly}+racon/log.txt",
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.prev_fa} {input.fq} > {output.sam} 2>{log}
        racon --threads {threads} --include-unpolished {input.fq} {output.sam} {input.prev_fa} > {output.fa} 2>>{log}
        ln -sr {output.fa} {output.link}
        """


# for running racon multiple iterations
rule raconX:
    conda:
        "env/conda-racon.yaml"
    threads: 5
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_{assembly}+racon{num}/output.fa",
        link="assemblies/{sample}_{assembly}+racon{num}.fa",
    log:
        "assemblies/{sample}_{assembly}+racon{num}/log.txt",
    shell:
        """
        DIR_temp=$(mktemp -d --suffix=.raconX)
        trap "rm -r $DIR_temp" EXIT

        cp {input.prev_fa} $DIR_temp/prev.fa

        for i in `seq 1 {wildcards.num}`
        do
            echo "Polishing round $i / {wildcards.num}" >>{log}
            minimap2 -ax map-ont -t {threads} $DIR_temp/prev.fa {input.fq} > $DIR_temp/map.sam 2>>{log}
            racon --threads {threads} --include-unpolished {input.fq} $DIR_temp/map.sam $DIR_temp/prev.fa > $DIR_temp/polished.fa 2>>{log}
            mv $DIR_temp/polished.fa $DIR_temp/prev.fa
        done

        mv $DIR_temp/prev.fa {output.fa}
        ln -sr {output.fa} {output.link}
        """


def get_model_for_sample(wildcards):
    if medaka_model is not None:
        return "-m " + medaka_model
    else:
        if map_medaka_model is not None:
            sample_base = wildcards.sample.split("+", 1)[0]
            if map_medaka_model.get(sample_base, False):
                return "-m " + map_medaka_model.get(sample_base, False)
            else:
                return ""
        else:
            return ""


rule medaka:
    conda:
        "env/conda-medaka.yaml"
    threads: 5
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq=get_ont_fq,
    output:
        fa="assemblies/{sample}_{assembly}+medaka/output.fa",
        link="assemblies/{sample}_{assembly}+medaka.fa",
    params:
        model=get_model_for_sample,
    log:
        "assemblies/{sample}_{assembly}+medaka/log.txt",
    shell:
        """
        medaka_consensus -f -i {input.fq} -d {input.prev_fa} -o assemblies/{wildcards.sample}_{wildcards.assembly}+medaka -t {threads} {params.model} >{log} 2>&1
        mv assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/consensus.fasta assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/output.fa
        ln -sr {output.fa} {output.link}
        """


rule pilon:
    conda:
        "env/conda-pilon.yaml"
    threads: 5
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq1=get_R1_fq,
        fq2=get_R2_fq,
    output:
        bam="assemblies/{sample}_{assembly}+pilon/map.bam",
        fa="assemblies/{sample}_{assembly}+pilon/output.fa",
        link="assemblies/{sample}_{assembly}+pilon.fa",
    log:
        "assemblies/{sample}_{assembly}+pilon/log.txt",
    shell:
        """
        bwa index -p assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.prev_fa} >{log} 2>&1
        bwa mem -t {threads} assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.fq1} {input.fq2} 2>>{log} | samtools sort -o {output.bam} - 2>>{log}
        samtools index {output.bam} 2>>{log}
        pilon -Xmx60G --genome {input.prev_fa} --frags {output.bam} --outdir assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/ --output pilon --changes --vcf >>{log} 2>&1
        mv assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/pilon.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


rule polca:
    conda:
        "env/conda-masurca.yaml"
    threads: 5
    shadow:
        "minimal"
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq1=get_R1_fq,
        fq2=get_R2_fq,
    output:
        fa="assemblies/{sample}_{assembly}+polca/output.fa",
        link="assemblies/{sample}_{assembly}+polca.fa",
    log:
        "assemblies/{sample}_{assembly}+polca/log.txt",
    shell:
        """
        polca.sh -t {threads} -a {input.prev_fa} -r '{input.fq1} {input.fq2}' >{log} 2>&1
        mv output.fa.PolcaCorrected.fa {output.fa}
        ln -sr {output.fa} {output.link}
        """


rule polypolish:
    conda:
        "env/conda-polypolish.yaml"
    threads: 5
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        fq1=get_R1_fq,
        fq2=get_R2_fq,
    output:
        fa="assemblies/{sample}_{assembly}+polypolish/output.fa",
        link="assemblies/{sample}_{assembly}+polypolish.fa",
    log:
        "assemblies/{sample}_{assembly}+polypolish/log.txt",
    shell:
        """
        bwa index -p assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/bwa_index {input.prev_fa} >{log} 2>&1
        bwa mem -a -t {threads} assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/bwa_index {input.fq1} > assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/alignments_R1.sam 2>>{log}
        bwa mem -a -t {threads} assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/bwa_index {input.fq2} > assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/alignments_R2.sam 2>>{log}
        polypolish_insert_filter.py --in1 assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/alignments_R1.sam --in2 assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/alignments_R2.sam --out1 assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/filtered_R1.sam --out2 assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/filtered_R2.sam >>{log} 2>&1
        polypolish {input.prev_fa} assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/filtered_R1.sam assemblies/{wildcards.sample}_{wildcards.assembly}+polypolish/filtered_R2.sam > {output.fa} 2>>{log}
        ln -sr {output.fa} {output.link}
        """


rule homopolish:
    conda:
        "env/conda-homopolish.yaml"
    threads: 1
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        ref="references/{ref}.fa",
    output:
        fa="assemblies/{sample}_{assembly}+homopolish/output_{ref}.fa",
        link="assemblies/{sample}_{assembly}+homopolish{ref}.fa",
    log:
        "assemblies/{sample}_{assembly}+homopolish/{ref}_log.txt",
    shell:
        """
        DIR_temp=$(mktemp -d --suffix=.raconX)
        trap "rm -r $DIR_temp" EXIT
        homopolish polish -a {input.prev_fa} -m R9.4.pkl -o $DIR_temp -l {input.ref} >{log} 2>&1
        cp $DIR_temp/*_homopolished.fasta {output.fa}
        ln -sr {output.fa} {output.link}
        """


rule proovframe_diamond_index:
    conda:
        "env/conda-proovframe.yaml"
    threads: 5
    input:
        "references-protein/{ref}.faa",
    output:
        "references-protein/{ref}.dmnd",
    log:
        "references-protein/{ref}-diamond-index.txt",
    shell:
        """
        diamond makedb -p {threads} --in {input} --db {output} >{log} 2>&1
        """


rule proovframe:
    conda:
        "env/conda-proovframe.yaml"
    threads: 1
    input:
        prev_fa="assemblies/{sample}_{assembly}/output.fa",
        ref="references-protein/{ref}.dmnd",
    output:
        fa="assemblies/{sample}_{assembly}+proovframe/output_{ref}.fa",
        tsv="assemblies/{sample}_{assembly}+proovframe/output_{ref}.tsv",
        link="assemblies/{sample}_{assembly}+proovframe{ref}.fa",
    log:
        "assemblies/{sample}_{assembly}+proovframe/{ref}_log.txt",
    shell:
        """
        proovframe map -d {input.ref} -o {output.tsv} {input.prev_fa} >{log} 2>&1
        proovframe fix -o {output.fa} {input.prev_fa} {output.tsv} >>{log} 2>&1
        ln -sr {output.fa} {output.link}
        """
