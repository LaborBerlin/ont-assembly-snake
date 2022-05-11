from glob import glob
import pandas as pd
shell.executable("/bin/bash")


filtlong_min_read_length = "1000"
if config.get('filtlong_min_read_length',False):
  filtlong_min_read_length = config['filtlong_min_read_length']
print("filtlong min. read length = " + filtlong_min_read_length)

#not needed anymore in flye 2.8
#genome_size = "6m"
#if config.get('genome_size',False):
#  genome_size = config['genome_size']
#print("genome size = " + genome_size)

medaka_model = None
if config.get('medaka_model',False):
  medaka_model = config['medaka_model']
  print("Medaka model = " + medaka_model)

map_medaka_model = None
if config.get('map_medaka_model',False):
  medaka_model_file = config['map_medaka_model']
  table = pd.read_csv(medaka_model_file, sep='\t', lineterminator='\n', header=None)
  map_medaka_model = dict(zip(table[0], table[1]))

wildcard_constraints:
  sample = "[^_]+",
  assembly = "[^_]+",
  sample_assembly = "[^/]+",
	num = "[0-9]+"

references, = glob_wildcards("references/{ref,[^/\\\\]+}.fa")

references_protein, = glob_wildcards("references-protein/{ref,[^/\\\\]+}.faa")

sample_assemblies, = glob_wildcards("assemblies/{sample_assembly,[^/]+}/")
#ignore symlinks in assemblies/folder, e.g. sample_flye.fa -> assemblies/sample_flye/output.fa
sample_assemblies = [a for a in sample_assemblies if not re.search('\.fa', a)]

# if any desired assembly requires homopolish then at least one reference genome should be provided
if not references and [string for string in sample_assemblies if "homopolish" in string]:
	quit("Error: must provide at least one reference genome sequence when using homopolish")

# if any desired assembly requires proovframe then at least one reference proteome should be provided
if not references_protein and [string for string in sample_assemblies if "proovframe" in string]:
	quit("Error: must provide at least one reference protein file when using proovframe")

list_outputs = expand("assemblies/{sample_assembly}/output.fa", sample_assembly = sample_assemblies)
list_outputs_links = expand("assemblies/{sample_assembly}.fa", sample_assembly = sample_assemblies)
# remove homopolish and proovframe assemblies from default list
list_outputs = [i for i in list_outputs if not re.search('homopolish|proovframe', i, re.IGNORECASE)]
list_outputs_links = [i for i in list_outputs_links if not re.search('homopolish|proovframe', i, re.IGNORECASE)]

# make lists for homopolish, one entry for each reference genome
list_outputs_homopolish = expand("assemblies/{sample_assembly}/output_{ref}.fa", ref = references, sample_assembly = [i for i in sample_assemblies if re.search('homopolish$', i, re.IGNORECASE)])
list_outputs_links_homopolish = expand("assemblies/{sample_assembly}{ref}.fa", ref = references, sample_assembly =  [i for i in sample_assemblies if re.search('homopolish$', i, re.IGNORECASE)])

# make lists for proovframe, one entry for each reference protein file
list_outputs_proovframe = expand("assemblies/{sample_assembly}/output_{ref}.fa", ref = references_protein, sample_assembly = [i for i in sample_assemblies if re.search('proovframe$', i, re.IGNORECASE)])
list_outputs_links_proovframe = expand("assemblies/{sample_assembly}{ref}.fa", ref = references_protein, sample_assembly =  [i for i in sample_assemblies if re.search('proovframe$', i, re.IGNORECASE)])

rule all:
	input:
		list_outputs,
		list_outputs_links,
		list_outputs_homopolish,
		list_outputs_links_homopolish,
		list_outputs_proovframe,
		list_outputs_links_proovframe

rule filtlong:
	threads: 1
	input:
		"fastq-ont/{sample}.fastq"
	output:
		"fastq-ont/{sample}+filtlong.fastq"
	log: "fastq-ont/{sample}_filtlong_log.txt"
	shell:
		"""
		filtlong --min_length {filtlong_min_read_length} {input} > {output} 2>{log}
		"""

rule filtlongX:
	threads: 1
	input:
		"fastq-ont/{sample}.fastq"
	output:
		"fastq-ont/{sample}+filtlong{num}.fastq"
	log: "fastq-ont/{sample}_filtlong{num}_log.txt"
	shell:
		"""
		filtlong --min_length {filtlong_min_read_length} -t {wildcards.num}000000 {input} > {output} 2>{log}
		"""

rule filtlongMql:
	threads: 1
	wildcard_constraints:
		mb = "[0-9]+",
		qweight = "[0-9]+",
		lweight = "[0-9]+"
	input:
		"fastq-ont/{sample}.fastq"
	output:
		"fastq-ont/{sample}+filtlong{mb},{qweight},{lweight}.fastq"
	log: "fastq-ont/{sample}_filtlong{mb},{qweight},{lweight}_log.txt"
	shell:
		"""
		filtlong --min_length {filtlong_min_read_length} --mean_q_weight {wildcards.qweight} --length_weight {wildcards.lweight}  -t {wildcards.mb}000000 {input} > {output} 2>{log}
		"""

rule filtlongMqln:
	threads: 1
	wildcard_constraints:
		mb = "[0-9]+",
		readlen = "[0-9]+",
		qweight = "[0-9]+",
		lweight = "[0-9]+"
	input:
		"fastq-ont/{sample}.fastq"
	output:
		"fastq-ont/{sample}+filtlong{mb},{qweight},{lweight},{readlen}.fastq"
	log: "fastq-ont/{sample}_filtlong{mb},{qweight},{lweight},{readlen}_log.txt"
	shell:
		"""
		filtlong --min_length {wildcards.readlen} --mean_q_weight {wildcards.qweight} --length_weight {wildcards.lweight}  -t {wildcards.mb}000000 {input} > {output} 2>{log}
		"""

rule unicycler:
  conda: "env/conda-unicycler.yaml"
	threads: 10
	input:
		fqont = "fastq-ont/{sample}.fastq",
		fq1 = "fastq-illumina/{sample}_R1.fastq",
		fq2 = "fastq-illumina/{sample}_R2.fastq"
	output:
		fa = "assemblies/{sample}_unicycler/output.fa",
		link = "assemblies/{sample}_unicycler.fa"
	log: "assemblies/{sample}_unicycler/log.txt"
	shell:
		"""
		unicycler -1 {input.fq1} -2 {input.fq2} -l {input.fqont} -t {threads} --keep 0 -o assemblies/{wildcards.sample}_unicycler/ >{log} 2>&1
		cp assemblies/{wildcards.sample}_unicycler/assembly.fasta {output.fa} 2>>{log}
		ln -sr {output.fa} {output.link}
		"""


#flye with default number of polishing rounds (=1 in flye v2.9)
rule flye:
  conda: "env/conda-flye.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_flye/output.fa",
		link = "assemblies/{sample}_flye.fa"
	log: "assemblies/{sample}_flye/log.txt"
	shell:
		"""
		flye --nano-raw {input.fq} -o assemblies/{wildcards.sample}_flye/ -t {threads} 2>{log}
		mv assemblies/{wildcards.sample}_flye/assembly.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

rule flyeX:
  conda: "env/conda-flye.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_flye{num}/output.fa",
		link = "assemblies/{sample}_flye{num}.fa"
	log: "assemblies/{sample}_flye{num}/log.txt"
	shell:
		"""
		flye --nano-raw {input.fq} -o assemblies/{wildcards.sample}_flye{wildcards.num}/ -t {threads} -i {wildcards.num} 2>{log}
		mv assemblies/{wildcards.sample}_flye{wildcards.num}/assembly.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

#flye for HQ ONT reads (from flye 2.9), with default number of polishing rounds (=1 in flye v2.9)
rule flyeHQ:
  conda: "env/conda-flye.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_flyeHQ/output.fa",
		link = "assemblies/{sample}_flyeHQ.fa"
	log: "assemblies/{sample}_flyeHQ/log.txt"
	shell:
		"""
		flye --nano-hq {input.fq} -o assemblies/{wildcards.sample}_flyeHQ/ -t {threads} 2>{log}
		mv assemblies/{wildcards.sample}_flye/assembly.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

rule flyeHQX:
  conda: "env/conda-flye.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_flyeHQ{num}/output.fa",
		link = "assemblies/{sample}_flyeHQ{num}.fa"
	log: "assemblies/{sample}_flyeHQ{num}/log.txt"
	shell:
		"""
		flye --nano-hq {input.fq} -o assemblies/{wildcards.sample}_flyeHQ{wildcards.num}/ -t {threads} -i {wildcards.num} 2>{log}
		mv assemblies/{wildcards.sample}_flyeHQ{wildcards.num}/assembly.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

#for running raven with default number of racon-polishing rounds (=2 in raven v0.0.8)
rule raven:
  conda: "env/conda-raven.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_raven/output.fa",
		link = "assemblies/{sample}_raven.fa"
	log: "assemblies/{sample}_raven/log.txt"
	shell:
		"""
		raven --disable-checkpoints -t {threads} {input.fq} >{output.fa} 2>{log}
		ln -sr {output.fa} {output.link}
		"""

#for running raven with racon polishing X times
rule ravenX:
  conda: "env/conda-raven.yaml"
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_raven{num}/output.fa",
		link = "assemblies/{sample}_raven{num}.fa"
	log: "assemblies/{sample}_raven{num}/log.txt"
	shell:
		"""
		raven --disable-checkpoints -p {wildcards.num} -t {threads} {input.fq} >{output.fa} 2>{log}
		ln -sr {output.fa} {output.link}
		"""

#for running racon once
rule racon:
  conda: "env/conda-racon.yaml"
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+racon/output.fa",
		link = "assemblies/{sample}_{assembly}+racon.fa",
		sam = temp("assemblies/{sample}_{assembly}+racon/map.sam")
	log: "assemblies/{sample}_{assembly}+racon/log.txt"
	shell:
		"""
		minimap2 -ax map-ont -t {threads} {input.prev_fa} {input.fq} > {output.sam} 2>{log}
		racon --threads {threads} --include-unpolished {input.fq} {output.sam} {input.prev_fa} > {output.fa} 2>>{log}
		ln -sr {output.fa} {output.link}
		"""

#for running racon multiple iterations
rule raconX:
  conda: "env/conda-racon.yaml"
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+racon{num}/output.fa",
		link = "assemblies/{sample}_{assembly}+racon{num}.fa"
	log: "assemblies/{sample}_{assembly}+racon{num}/log.txt"
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
      sample_base = wildcards.sample.split("+",1)[0]
      if map_medaka_model.get(sample_base,False):
        return "-m " + map_medaka_model.get(sample_base,False)
      else:
        return ""
    else:
      return ""

rule medaka:
  conda: "env/conda-medaka.yaml"
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+medaka/output.fa",
		link = "assemblies/{sample}_{assembly}+medaka.fa"
	params:
		model = get_model_for_sample
	log: "assemblies/{sample}_{assembly}+medaka/log.txt"
	message: "Medaka: {wildcards.sample}, {wildcards.assembly}, {params.model}"
	shell:
		"""
		medaka_consensus -f -i {input.fq} -d {input.prev_fa} -o assemblies/{wildcards.sample}_{wildcards.assembly}+medaka -t {threads} {params.model} >{log} 2>&1
		mv assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/consensus.fasta assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/output.fa
		ln -sr {output.fa} {output.link}
		"""

rule pilon:
  conda: "env/conda-pilon.yaml"
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq1 = "fastq-illumina/{sample}_R1.fastq",
		fq2 = "fastq-illumina/{sample}_R2.fastq"
	output:
		bam = "assemblies/{sample}_{assembly}+pilon/map.bam",
		fa = "assemblies/{sample}_{assembly}+pilon/output.fa",
		link = "assemblies/{sample}_{assembly}+pilon.fa"
	log: "assemblies/{sample}_{assembly}+pilon/log.txt"
	shell:
		"""
		bwa index -p assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.prev_fa} >{log} 2>&1
		bwa mem -t {threads} assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.fq1} {input.fq2} 2>>{log} | samtools sort -o {output.bam} - 2>>{log}
		samtools index {output.bam} 2>>{log}
		pilon -Xmx60G --genome {input.prev_fa} --frags {output.bam} --outdir assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/ --output pilon --changes --vcf >>{log} 2>&1
		mv assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/pilon.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

rule polypolish:
  conda: "env/conda-polypolish.yaml"
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq1 = "fastq-illumina/{sample}_R1.fastq",
		fq2 = "fastq-illumina/{sample}_R2.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+polypolish/output.fa",
		link = "assemblies/{sample}_{assembly}+polypolish.fa"
	log: "assemblies/{sample}_{assembly}+polypolish/log.txt"
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
  conda: "env/conda-homopolish.yaml"
	threads: 1
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		ref = "references/{ref}.fa"
	output:
		fa = "assemblies/{sample}_{assembly}+homopolish/output_{ref}.fa",
		link = "assemblies/{sample}_{assembly}+homopolish{ref}.fa"
	log: "assemblies/{sample}_{assembly}+homopolish/{ref}_log.txt"
	shell:
		"""
		DIR_temp=$(mktemp -d --suffix=.raconX)
		trap "rm -r $DIR_temp" EXIT
		homopolish polish -a {input.prev_fa} -m R9.4.pkl -o $DIR_temp -l {input.ref} >{log} 2>&1
		cp $DIR_temp/*_homopolished.fasta {output.fa}
		ln -sr {output.fa} {output.link}
		"""

rule proovframe_diamond_index:
  conda: "env/conda-proovframe.yaml"
	threads: 5
	input: "references-protein/{ref}.faa"
	output: "references-protein/{ref}.dmnd"
	log: "references-protein/{ref}-diamond-index.txt"
	shell:
		"""
		diamond makedb -p {threads} --in {input} --db {output} >{log} 2>&1
		"""

rule proovframe:
  conda: "env/conda-proovframe.yaml"
	threads: 1
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		ref = "references-protein/{ref}.dmnd"
	output:
		fa = "assemblies/{sample}_{assembly}+proovframe/output_{ref}.fa",
		tsv = "assemblies/{sample}_{assembly}+proovframe/output_{ref}.tsv",
		link = "assemblies/{sample}_{assembly}+proovframe{ref}.fa"
	log: "assemblies/{sample}_{assembly}+proovframe/{ref}_log.txt"
	shell:
		"""
		proovframe map -d {input.ref} -o {output.tsv} {input.prev_fa} >{log} 2>&1
		proovframe fix -o {output.fa} {input.prev_fa} {output.tsv} >{log} 2>&1
		ln -sr {output.fa} {output.link}
		"""

