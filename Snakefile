from glob import glob
shell.executable("/bin/bash")

genome_size = "6m"
if config.get('genome_size',False):
  genome_size = config['genome_size']

print(genome_size)

medaka_model = "r941_min_high_g344"
if config.get('medaka_model',False):
  medaka_model = config['medaka_model']

flye_iterations = 4
if config.get('flye_iterations',False):
  flye_iterations = config['flye_iterations']

wildcard_constraints:
  sample = "[^_]+",
  assembly = "[^_]+",
  sample_assembly = "[^/]+",
	num = "[1-9]"

sample_assemblies, = glob_wildcards("assemblies/{sample_assembly,[^/]+}/")

list_outputs = expand("assemblies/{sample_assembly}/output.fa", sample_assembly = sample_assemblies)

rule all:
	input:
		list_outputs

rule flye:
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq",
	output:
		fa = "assemblies/{sample}_flye/output.fa"
	log: "assemblies/{sample}_flye/log.txt"
	shell:
		"""
		flye --nano-raw {input.fq} -g {genome_size} -o assemblies/{wildcards.sample}_flye/ -t {threads} -i {flye_iterations} 2>{log}
		mv assemblies/{wildcards.sample}_flye/assembly.fasta {output.fa}
		"""

rule raven:
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_raven/output.fa"
	log: "assemblies/{sample}_raven/log.txt"
	shell:
		"""
		raven -t {threads} {input.fq} >{output.fa} 2>{log}
		"""

#for running raven with racon polishing
rule ravenX:
	threads: 5
	input:
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_raven{num}/output.fa"
	log: "assemblies/{sample}_raven{num}/log.txt"
	shell:
		"""
		raven -p {wildcards.num} -t {threads} {input.fq} >{output.fa} 2>{log}
		"""

#for running racon once
rule racon:
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+racon/output.fa",
		sam = temp("assemblies/{sample}_{assembly}+racon/map.sam")
	log: "assemblies/{sample}_{assembly}+racon/log.txt"
	shell:
		"""
		minimap2 -ax map-ont -t {threads} {input.prev_fa} {input.fq} > {output.sam} 2>{log}
		racon --threads {threads} --include-unpolished {input.fq} {output.sam} {input.prev_fa} > {output.fa}
		"""

#for running racon multiple iterations
rule raconX:
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		fa = "assemblies/{sample}_{assembly}+racon{num}/output.fa"
	log: "assemblies/{sample}_{assembly}+racon{num}/log.txt"
	shell:
		"""
		DIR_temp=$(mktemp -d --suffix=.raconnn)
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
		"""

rule medaka:
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq = "fastq-ont/{sample}.fastq"
	output:
		"assemblies/{sample}_{assembly}+medaka/output.fa"
	log: "assemblies/{sample}_{assembly}+medaka/log.txt"
	shell:
		"""
		medaka_consensus -i {input.fq} -d {input.prev_fa} -o assemblies/{wildcards.sample}_{wildcards.assembly}+medaka -t {threads} -m {medaka_model} >{log} 2>&1
		mv assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/consensus.fasta assemblies/{wildcards.sample}_{wildcards.assembly}+medaka/output.fa
		"""

rule pilon:
	threads: 5
	input:
		prev_fa = "assemblies/{sample}_{assembly}/output.fa",
		fq1 = "fastq-illumina/{sample}_R1.fastq",
		fq2 = "fastq-illumina/{sample}_R2.fastq"
	output:
		bam = "assemblies/{sample}_{assembly}+pilon/map.bam",
		fa = "assemblies/{sample}_{assembly}+pilon/output.fa"
	log: "assemblies/{sample}_{assembly}+pilon/log.txt"
	shell:
		"""
		bwa index -p assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.prev_fa} >{log} 2>&1
		bwa mem -t {threads} assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/bwa_index {input.fq1} {input.fq2} 2>>{log} | samtools sort -o {output.bam} - 2>>{log}
		samtools index {output.bam} 2>>{log}
		pilon -Xmx60G --genome {input.prev_fa} --frags {output.bam} --outdir assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/ --output pilon --changes --vcf >>{log} 2>&1
		mv assemblies/{wildcards.sample}_{wildcards.assembly}+pilon/pilon.fasta {output.fa}
		"""

