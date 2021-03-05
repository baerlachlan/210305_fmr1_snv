#!/bin/bash

## This script follows the recommended gatk workflow for RNA-seq short variant discovery (SNPs + Indels):
## - https://gatk.broadinstitute.org/hc/en-us/articles/360035531192-RNAseq-short-variant-discovery-SNPs-Indels-

## Minor adjustments have been made as follows:
## - The workflow assumes raw data is in the unmapped BAM format
##   - As the raw data for this dataset is in FASTQ format, the initial SamToFastq method was skipped
##   - Later in the workflow, the MergeBamAlignment step was also skipped
##     - Unmapped BAM files contain metadata and this step allows for merging of that metadata into the aligned BAMS
##     - This was not needed as raw data was FASTQ which does not contain metadata
## - A bootstrapping approach was taken to define a set of known variants
##   - Our zebrafish are most likely quite different to the reference models due to years of breeding
##   - It is better to therefore generate a set of known variants from our data than to use pre-defined ones in a database (dbSNP , Ensembl)
##   - The reasoning and process for this is located at https://gatk.broadinstitute.org/hc/en-us/articles/360035890531?id=44 (Section 3)

SAMPLES = \
["A", "D", "G", "L", "S2", "S4", "S5", "S8"]

refs = "/hpcfs/users/a1647910/refs/ensembl-release-94/danio_rerio/"
readLength = 150

rule all:
	input:
		# refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict",
		# refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
		expand("00_rawData/FastQC/{SAMPLE}_1_fastqc.zip", SAMPLE = SAMPLES),
		expand("00_rawData/FastQC/{SAMPLE}_1_fastqc.html", SAMPLE = SAMPLES),
		expand("00_rawData/FastQC/{SAMPLE}_2_fastqc.zip", SAMPLE = SAMPLES),
		expand("00_rawData/FastQC/{SAMPLE}_2_fastqc.html", SAMPLE = SAMPLES),
		expand("01_trimmedData/FastQC/{SAMPLE}_1_fastqc.zip", SAMPLE = SAMPLES),
		expand("01_trimmedData/FastQC/{SAMPLE}_1_fastqc.html", SAMPLE = SAMPLES),
		expand("01_trimmedData/FastQC/{SAMPLE}_2_fastqc.zip", SAMPLE = SAMPLES),
		expand("01_trimmedData/FastQC/{SAMPLE}_2_fastqc.html", SAMPLE = SAMPLES),
		# expand("02_alignedData/FastQC/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.zip", SAMPLE = SAMPLES),
		# expand("02_alignedData/FastQC/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.html", SAMPLE = SAMPLES),
		# expand("03_markDuplicates/FastQC/{SAMPLE}_fastqc.zip", SAMPLE = SAMPLES),
		# expand("03_markDuplicates/FastQC/{SAMPLE}_fastqc.html", SAMPLE = SAMPLES),
		# expand("11_filterVariants/vcf/{SAMPLE}.vcf.gz", SAMPLE = SAMPLES),
		# expand("11_filterVariants/vcf/{SAMPLE}.vcf.gz.tbi", SAMPLE = SAMPLES)
		bam = "03_markDuplicates/bam/A.bam",
		bamIndex = "03_markDuplicates/bam/A.bai",
		metrics = "03_markDuplicates/log/A.metrics"


rule fastqc_raw:
	input:
		"00_rawData/fastq/{SAMPLE}.fq.gz"
	output:
		"00_rawData/FastQC/{SAMPLE}_fastqc.zip",
		"00_rawData/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "00_rawData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule trim:
	input:
		R1 = "00_rawData/fastq/{SAMPLE}_1.fq.gz",
		R2 = "00_rawData/fastq/{SAMPLE}_2.fq.gz"
	output:
		R1 = "01_trimmedData/fastq/{SAMPLE}_1.fq.gz",
		R2 = "01_trimmedData/fastq/{SAMPLE}_2.fq.gz",
		setting = "01_trimmedData/fastq/{SAMPLE}.settings",
		discard = "01_trimmedData/fastq/{SAMPLE}.discarded.gz",
		single = "01_trimmedData/fastq/{SAMPLE}.singleton.truncated.gz" ## Only for paired end
	params:
		bname = "01_trimmedData/fastq/{SAMPLE}"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 2,
		mem_mb = 2000,
		hours = 2,
		mins = 0
	shell:
		"""
		AdapterRemoval --gzip \
		--file1 {input.R1} \
		--file2 {input.R2} \
		--output1 {output.R1} \
		--output2 {output.R2} \
		--basename {params.bname} \
		--trimns \
		--trimqualities \
		--minquality 20 \
		--minlength 35 \
		--threads {resources.cpu}
		"""

rule fastqc_trim:
	input:
		"01_trimmedData/fastq/{SAMPLE}.fq.gz"
	output:
		"01_trimmedData/FastQC/{SAMPLE}_fastqc.zip",
		"01_trimmedData/FastQC/{SAMPLE}_fastqc.html"
	params:
		outDir = "01_trimmedData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule star_index:
	input:
		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
		gtf = refs + "Danio_rerio.GRCz11.94.chr.gtf.gz"
	output:
		directory(refs + "star/")
	params:
		overhang = readLength-1
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		hours = 0,
		mins = 30
	shell:
		"""
		zcat {input.gtf} > temp.gtf

		STAR \
		--runThreadN {resources.cpu} \
		--runMode genomeGenerate \
		--genomeDir {output} \
		--genomeFastaFiles {input.refFa} \
		--sjdbGTFfile temp.gtf \
		--sjdbOverhang {params.overhang}

		rm temp.gtf
		"""

rule align:
	input:
		R1 = "01_trimmedData/fastq/{SAMPLE}_1.fq.gz",
		R2 = "01_trimmedData/fastq/{SAMPLE}_2.fq.gz",
		starIndex = refs + "star/"
	output:
		"02_alignedData/bam/{SAMPLE}_Aligned.sortedByCoord.out.bam"
	params:
		overhang = readLength-1,
		bname = "02_alignedData/bam/{SAMPLE}_"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 16,
		ntasks = 1,
		mem_mb = 32000,
		hours = 5,
		mins = 0
	shell:
		"""
		STAR \
		--genomeDir {input.starIndex}\
		--runThreadN {resources.cpu} \
		--readFilesIn {input.R1} {input.R2} \
		--readFilesCommand "gunzip -c" \
		--sjdbOverhang {params.overhang} \
		--outSAMtype BAM SortedByCoordinate \
		--twopassMode Basic \
		--outFileNamePrefix {params.bname}


		mkdir -p 02_alignedData/log
		mv {params.bname}*out 02_alignedData/log
		mv {params.bname}*tab 02_alignedData/log
		"""

rule fastqc_align:
	input:
		"02_alignedData/bam/{SAMPLE}_Aligned.sortedByCoord.out.bam"
	output:
		"02_alignedData/FastQC/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.zip",
		"02_alignedData/FastQC/{SAMPLE}_Aligned.sortedByCoord.out_fastqc.html"
	params:
		outDir = "02_alignedData/FastQC/"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 1,
		ntasks = 1,
		mem_mb = 2000,
		hours = 1,
		mins = 0
	shell:
		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

rule mark_duplicates:
	input:
		# "02_alignedData/bam/{SAMPLE}_Aligned.sortedByCoord.out.bam"
		"02_alignedData/bam/A_Aligned.sortedByCoord.out.bam"
	output:
		# bam = "03_markDuplicates/bam/{SAMPLE}.bam",
		# bamIndex = "03_markDuplicates/bam/{SAMPLE}.bai",
		# metrics = "03_markDuplicates/log/{SAMPLE}.metrics"
		bam = "03_markDuplicates/bam/A.bam",
		bamIndex = "03_markDuplicates/bam/A.bai",
		metrics = "03_markDuplicates/log/A.metrics"
	conda:
		"snakemake/envs/default.yaml"
	resources:
		cpu = 2,
		ntasks = 1,
		mem_mb = 4000,
		hours = 4,
		mins = 0
	shell:
		"""
		gatk \
 	    	MarkDuplicates \
 	        --INPUT {input} \
 	        --OUTPUT {output.bam}  \
 	        --CREATE_INDEX true \
 	        --VALIDATION_STRINGENCY SILENT \
 	        --METRICS_FILE {output.metrics}
		"""

# rule fastqc_duplicates:
# 	input:
# 		"03_markDuplicates/bam/{SAMPLE}.bam"
# 	output:
# 		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.zip",
# 		"03_markDuplicates/FastQC/{SAMPLE}_fastqc.html"
# 	params:
# 		outDir = "03_markDuplicates/FastQC/"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 1,
# 		mem_mb = 2000,
# 		hours = 1,
# 		mins = 0
# 	shell:
# 		"fastqc -t {resources.cpu} -o {params.outDir} --noextract {input}"

# ## Reference dictionary and index needs to be created as described in:
# ## https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format
# rule ref_dict:
# 	input:
# 		refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
# 	output:
# 		refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 1,
# 		mem_mb = 1000,
# 		hours = 0,
# 		mins = 10
# 	shell:
# 		"gatk CreateSequenceDictionary -R {input}"

# rule ref_index:
# 	input:
# 		refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa"
# 	output:
# 		refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 1,
# 		mem_mb = 200,
# 		hours = 0,
# 		mins = 10
# 	shell:
# 		"samtools faidx {input}"

# rule splitNCigar:
# 	input:
# 		bam = "03_markDuplicates/bam/{SAMPLE}.bam",
# 		bamIndex = "03_markDuplicates/bam/{SAMPLE}.bai",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	output:
# 		bam = "04_splitNCigar/bam/{SAMPLE}.bam",
# 		bamIndex = "04_splitNCigar/bam/{SAMPLE}.bai"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 8,
# 		ntasks = 1,
# 		mem_mb = 16000,
# 		hours = 10,
# 		mins = 0
# 	shell:
# 		"""
# 		gatk \
#             SplitNCigarReads \
#             -R {input.refFa} \
#             -I {input.bam} \
#             -O {output.bam}
# 		"""

# ## gatk HaplotypeCaller requires read group (RG) tags
# ## An explanation of this is found at https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
# rule addRG:
# 	input:
# 		bam = "04_splitNCigar/bam/{SAMPLE}.bam",
# 		bamIndex = "04_splitNCigar/bam/{SAMPLE}.bai"
# 	output:
# 		bam = "05_addRG/bam/{SAMPLE}.bam",
# 		bamIndex = "05_addRG/bam/{SAMPLE}.bai"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 1000,
# 		hours = 0,
# 		mins = 30
# 	shell:
# 		"""
# 		gatk \
# 			AddOrReplaceReadGroups \
#     		-I {input} \
#    			-O {output.bam} \
#     		-SORT_ORDER coordinate \
#     		-RGID 1 \
#     		-RGLB SORL1 \
#     		-RGSM {SAMPLE} \
# 			-RGPU SORL1 \
#     		-RGPL illumina \
#     		-CREATE_INDEX True
# 		"""

# ## This step allows for defining a known set of variants based on the data, opposed to using a pre-existing database
# rule callVariants_noRecal:
# 	input:
# 		bam = "05_addRG/bam/{SAMPLE}.bam",
# 		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	output:
# 		vcf = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz.tbi"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 6500,
# 		hours = 5,
# 		mins = 0
# 	shell:
# 		"""
# 		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
# 		HaplotypeCaller \
# 		-R {input.refFa} \
# 		-I {input.bam} \
# 		-O {output.vcf} \
# 		-dont-use-soft-clipped-bases \
# 		--standard-min-confidence-threshold-for-calling 20
# 		"""

# rule knownVariants:
# 	input:
# 		vcf = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "06_callVariants_noRecal/vcf/{SAMPLE}.vcf.gz.tbi",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	output:
# 		vcf = "07_knownVariants/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "07_knownVariants/vcf/{SAMPLE}.vcf.gz.tbi"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 1000,
# 		hours = 0,
# 		mins = 10
# 	shell:
# 		"""
# 		gatk \
# 		    VariantFiltration \
# 			--R {input.refFa} \
# 			--V {input.vcf} \
# 			--window 35 \
# 			--cluster 3 \
# 			--filter-name "FS" \
# 			--filter "FS > 30.0" \
# 			--filter-name "QD" \
# 			--filter "QD < 2.0" \
# 			-O {output.vcf}
# 		"""

# ## We now follow the GATK workflow as normal
# rule baseRecal:
# 	input:
# 		bam = "05_addRG/bam/{SAMPLE}.bam",
# 		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict",
# 		dbsnp = "07_knownVariants/vcf/{SAMPLE}.vcf.gz",
# 		dbsnpIndex = "07_knownVariants/vcf/{SAMPLE}.vcf.gz.tbi"
# 	output:
# 		"08_baseRecal/{SAMPLE}_recal"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 4000,
# 		hours = 0,
# 		mins = 10
# 	shell:
# 		"""
# 		gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
#             -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
#             -Xloggc:gc_log.log -Xms4000m" \
#             BaseRecalibrator \
#             -R {input.refFa} \
#             -I {input.bam} \
#             --use-original-qualities \
#             -O {output} \
#             -known-sites {input.dbsnp}
# 		"""

# rule applyRecal:
# 	input:
# 		bam = "05_addRG/bam/{SAMPLE}.bam",
# 		bamIndex = "05_addRG/bam/{SAMPLE}.bai",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict",
# 		recal = "08_baseRecal/{SAMPLE}_recal"
# 	output:
# 		bam = "09_applyRecal/bam/{SAMPLE}.bam",
# 		bamIndex = "09_applyRecal/bam/{SAMPLE}.bai"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 4000,
# 		hours = 0,
# 		mins = 30
# 	shell:
# 		"""
# 		gatk \
#             --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
#             -XX:+PrintGCDetails -Xloggc:gc_log.log \
#             -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m" \
#             ApplyBQSR \
#             --add-output-sam-program-record \
#             -R {input.refFa} \
#             -I {input.bam} \
#             --use-original-qualities \
#             -O {output.bam} \
#             --bqsr-recal-file {input.recal}
# 		"""

# rule callVariants:
# 	input:
# 		bam = "09_applyRecal/bam/{SAMPLE}.bam",
# 		bamIndex = "09_applyRecal/bam/{SAMPLE}.bai",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	output:
# 		vcf = "10_callVariants/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "10_callVariants/vcf/{SAMPLE}.vcf.gz.tbi"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 6500,
# 		hours = 5,
# 		mins = 0
# 	shell:
# 		"""
# 		gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
# 		HaplotypeCaller \
# 		-R {input.refFa} \
# 		-I {input.bam} \
# 		-O {output.vcf} \
# 		-dont-use-soft-clipped-bases \
# 		--standard-min-confidence-threshold-for-calling 20
# 		"""

# rule filterVariants:
# 	input:
# 		vcf = "10_callVariants/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "10_callVariants/vcf/{SAMPLE}.vcf.gz.tbi",
# 		refFa = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa",
# 		refIndex = refs + "Danio_rerio.GRCz11.dna.primary_assembly.fa.fai",
# 		refDict = refs + "Danio_rerio.GRCz11.dna.primary_assembly.dict"
# 	output:
# 		vcf = "11_filterVariants/vcf/{SAMPLE}.vcf.gz",
# 		vcfIndex = "11_filterVariants/vcf/{SAMPLE}.vcf.gz.tbi"
# 	conda:
# 		"snakemake/envs/default.yaml"
# 	resources:
# 		cpu = 1,
# 		ntasks = 2,
# 		mem_mb = 1000,
# 		hours = 0,
# 		mins = 10
# 	shell:
# 		"""
# 		gatk \
# 		    VariantFiltration \
# 			--R {input.refFa} \
# 			--V {input.vcf} \
# 			--window 35 \
# 			--cluster 3 \
# 			--filter-name "FS" \
# 			--filter "FS > 30.0" \
# 			--filter-name "QD" \
# 			--filter "QD < 2.0" \
# 			-O {output.vcf}
# 		"""