import os
import pandas as pd


# samples = ["father", "mother", "proband"]
samples = ["proband", "mother", "father"]
aligners = ["bwa"]
variant_callers = ["GATK"]

RS = "R"
EXT = "fastq"

samples_dir = "/home/marc/WES_GATK/test/samples"
out_dir = "/home/marc/WES_GATK/test/snake_out2"
ref_bwa = "/home/marc/Desktop/data/refs/indexes/bwa/hg38/hg38"
ref_bwa_path = "/home/marc/Desktop/data/refs/indexes/bwa"
ref_fasta = "~/Desktop/data/refs/fa/hg38/hg38.fa"
ref_fasta_path = "/home/marc/Desktop/data/refs/fa/hg38"
bed_file = "/home/marc/Desktop/data/refs/BED/Twist_Exome_Core_Covered_Targets_hg38.bed"
known_variants = "/home/marc/Desktop/data/refs/vcf/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz"
nirvana_path = "/home/marc/Nirvana"
Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll"
Nirvana_supplementray = f"{nirvana_path}/Data/SupplementaryAnnotation/GRCh38/"
Nirvana_ref = f"{nirvana_path}/Data/References/Homo_sapiens.GRCh38.Nirvana.dat"
Nirvana_cache = f"{nirvana_path}/Data/Cache/GRCh38/Both"
gff = "/home/marc/Desktop/data/refs/gtf/Homo_sapiens.GRCh38.109.gff3.gz"
env_path = "/home/marc/WES_GATK/test_workflow/envs"

workdir: out_dir

def get_gvcf(wildcards):
    return [f"04_calling/{sample}_raw.gvcf.gz" \
        for sample in samples]

def get_bams(wildcards):
    return [f"03_bamPrep/{sample}.pqsr.bam" \
        for sample in samples]

def get_final_output(wildcards):
    final_output = []
    final_output.extend(expand(
        f"00_QC/{{sample}}_R{{R}}_fastqc.zip",
            sample = samples,
            R = [1,2]
    ))
    final_output.extend(expand(
        "02_alignment/QC/{sample}.cov",
            sample = samples
    ))
    final_output.extend(expand(
        "03_bamPrep/QC/{sample}_Qualimap",
        sample = samples
    ))
    final_output.extend(expand(
        "03_bamPrep/QC/{sample}.pdf",
            sample = samples
    ))
    final_output.extend(
        [
            "04_calling/QC/bcftools.stats",
            "04_calling/QC/bcftools_plots",
            "04_calling/QC/bcftools_csq.vcf",
            "04_calling/snvs/variants_genotyped.gvcf.gz",
            "04_calling/indels/variants_genotyped.gvcf.gz",
            "05_Annotation/Nirvana/Annotation.json.gz",
            "05_Annotation/ANNOVAR"
        ]
    )
    final_output.extend(expand(
        "04_calling/QC/{sample}.eval.grp",
            sample = samples
    ))



    return final_output

rule all:
    input:
        "multiqc/multiqc_report.html"


######################################
#####            QC         ##########
######################################


rule trimmomatic:
    input:
        R1 = f"{samples_dir}/{{sample}}_R1.fastq",
        R2 = f"{samples_dir}/{{sample}}_R2.fastq"

    
    conda: "GUAP"

    output:
        log="logs/trimmomatic/{sample}.log",
        summary="logs/trimmomatic/{sample}.summary",
        nf1 = f"01_trimmomatic/{{sample}}_R1.{EXT}",
        nf2 = f"01_trimmomatic/{{sample}}_R2.{EXT}",
        nfu1=temp(f"01_trimmomatic/U/{{sample}}_R1_U.{EXT}"),
        nfu2=temp(f"01_trimmomatic/U/{{sample}}_R2_U.{EXT}")

    benchmark: "benchamrks/QC/{sample}_trim.txt"
    threads: 4
    params:
        size = 4,
        quality = 10,
        extra = '',
        minlen = 50

    resources:
        mem_mb=5120,
        cores=4,
        mem_gb=5,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    log: 
        "logs/trimmomatic/{sample}.txt"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
                -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
                SLIDINGWINDOW:{params.size}:{params.quality} MINLEN:{params.minlen} > {log} 2>&1
        """

rule Fastqc:
    input:
        f"01_trimmomatic/{{sample}}_R{{R}}.{EXT}"
    log:
        f"logs/QC/QC_{{sample}}_R{{R}}.log"
    
    conda: "GUAP"

    output:
        zip=f"00_QC/{{sample}}_R{{R}}_fastqc.zip",
        html=f"00_QC/{{sample}}_R{{R}}_fastqc.html"
        
    benchmark: f"benchamrks/QC/{{sample}}_R{{R}}fastqc.txt"
    threads: 2
    params:
        path="00_QC"
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        "fastqc {input} --threads {threads} -o {params.path} > {log} 2>&1"

######################################
#####           align       ##########
######################################

rule index_ref:
    input: f"{ref_fasta}"
    
    conda: "GUAP"

    output: f"{ref_fasta}.fai"
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: "samtools faidx {input}"


rule refrence_dict:
    input: f"{ref_fasta_path}/hg38.fa"
    
    conda: "GUAP"

    output: f"{ref_fasta_path}/hg38.dict"
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: "picard CreateSequenceDictionary -R {input}"

rule bwa_index:
    input: f"{ref_fasta_path}/hg38.fa"
    
    conda: "GUAP"

    output: directory(f"{ref_bwa_path}/hg38")
    params:
        prefix = "hg38"
    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: 
        "bwa index -p {output}/{params.prefix} {input}"


rule bwa_align:
    input:
        R1 = f"01_trimmomatic/{{sample}}_R1.{EXT}",
        R2 = f"01_trimmomatic/{{sample}}_R2.{EXT}",

    
    conda: "GUAP"

    output:
        temp("02_alignment/{sample}.sam")

    threads: 4
    params:
        index = ref_bwa,
        fa = ref_fasta

    log: 
        bwa = "logs/{sample}_bwa.log",

    benchmark: "benchamrks/{sample}_bwa.txt"
    resources:
        mem_mb=32768,
        cores=4,
        mem_gb=32,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        R1={input.R1}
        SM=$(basename $R1 | cut -d"_" -f1)   
        PL="Illumina"
        LB="001"
        name=$(basename $R1 | cut -d'_' -f1)
        RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
        PU=$RGID.$LB 
        bwa mem -t {threads} -M \
            -R "@RG\\tID:$RGID\\tSM:$SM\\tPL:$PL\\tLB:$LB\\tPU:$PU" {params.index} {input.R1} {input.R2} > {output} 2> {log.bwa}
        """


rule sort_and_convert_sam:
    input:
        "02_alignment/{sample}.sam"

    
    conda: "GUAP"

    output:
        "02_alignment/{sample}.sorted.bam"
    
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        samtools sort {input} -o {output}
        samtools index {output}
        """


rule QC_alignment:
    input:
        "02_alignment/{sample}.sorted.bam"

    
    conda: "GUAP"

    output:
        cov = "02_alignment/QC/{sample}.cov",
        stats = "02_alignment/QC/{sample}.stats"
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > {output.cov}
        samtools flagstat {input} > {output.stats}
        """


######################################
#####      bam processing   ##########
######################################

rule mrk_duplicates:
    input:
        "02_alignment/{sample}.sorted.bam"

    
    conda: "GUAP"

    output:
        bam = "02_alignment/{sample}.dedub.bam",
        matrix = "02_alignment/{sample}.dedub.matrix"

    resources:
        mem_mb=4096,
        cores=4,
        mem_gb=4,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    log: 
        "logs/{sample}.dedub.log"

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} > {log}
        """

rule BaseRecalibrator:
    input: "02_alignment/{sample}.dedub.bam"
    
    conda: "GUAP"

    output: "03_bamPrep/{sample}.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta,

    threads: 4
    benchmark: "benchamrks/{sample}_GATK_pqsr.txt"
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """

rule applyBaseRecalibrator:
    input: 
        bam = "02_alignment/{sample}.dedub.bam",
        report = "03_bamPrep/{sample}.report"
    benchmark: "benchamrks/{sample}_GATK_apply_BQSR.txt"
    
    conda: "GUAP"

    output: "03_bamPrep/{sample}.pqsr.bam"
    threads: 1
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" ApplyBQSR  -R {params.ref} \
            -I {input.bam} --emit-original-quals \
            -bqsr {input.report} -O {output} \
            --add-output-sam-program-record
        """

rule bqsr_calibrated_report:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "GUAP"

    output: "03_bamPrep/{sample}_pqsr.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    threads: 1
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """

rule AnalyzeCovariates:
    input: 
        raw = "03_bamPrep/{sample}.report", 
        bqsr = "03_bamPrep/{sample}_pqsr.report"

    
    conda: "GUAP"

    output: "03_bamPrep/QC/{sample}.pdf"

    threads: 2
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" AnalyzeCovariates -before {input.raw} \
            -after {input.bqsr} -plots {output}
        """

rule qualimap:
    input:
        "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "GUAP"

    output:
        directory("03_bamPrep/QC/{sample}_Qualimap")

    threads: 2
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt  

    shell:
        """
        qualimap \
            bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            --genome-gc-distr HUMAN \
            -nt {threads} \
            -skip-duplicated \
            --skip-dup-mode 0 \
            -outdir {output} \
            -outformat HTML
        """


######################################
#####     variant calling    #########
######################################

rule HaplotypeCaller:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "GUAP"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file

    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" HaplotypeCaller -R {params.ref} \
            -L {params.bed} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output}
        """

rule combine_gvcf:
    input:
        gvcfs=get_gvcf
    
    conda: "GUAP"

    output:
        "04_calling/variants.gvcf.gz"
    params:
        gvcfs = lambda wildcards, input: [f" --variant {v}" for v in input["gvcfs"]],
        ref = ref_fasta

    threads: 1
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
            -O {output}
        """

rule genotype_gvcfs:
    input:
        "04_calling/variants.gvcf.gz"
    
    conda: "GUAP"

    output:
        "04_calling/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    threads: 4
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" GenotypeGVCFs \
            -R {params.ref} -V {input} -O {output}
        """

rule VariantEval:
    input: "04_calling/{sample}_raw.gvcf.gz"
    
    conda: "GUAP"

    output: "04_calling/QC/{sample}.eval.grp"
    threads:1
    params:
        ref = ref_fasta
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        touch {output}
        gatk VariantEval \
            -R {params.ref} \
            --eval:{wildcards.sample} {input} \
            -O {output}
        """

rule bcftools_stats:
    input: "04_calling/QC/bcftools_csq.vcf"
    
    conda: "wes"

    output: "04_calling/QC/bcftools.stats"
    threads: 1
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
            bcftools stats {input} > {output}
        """

rule plot_bcftools_stats:
    input: "04_calling/QC/bcftools.stats"
    
    conda: "wes"

    output: directory("04_calling/QC/bcftools_plots")
    threads: 1 
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        plot-vcfstats -p {output} {input}
        """

rule consequence:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output: "04_calling/QC/bcftools_csq.vcf"
    threads: 1
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads:1
    params:
        ref = ref_fasta,
        gff = gff

    shell:
        """
        bcftools csq -f {params.ref} -g {params.gff} {input} -Ov -o {output}
        """

######################################
#####   vairant annotation   #########
######################################


rule Nirvana:
    input:
        "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output:
        "05_Annotation/Nirvana/Annotation.json.gz"
    threads: 1
    resources:
        mem_mb=8192,
        cores=1,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        Nirvana_supplementray = f"{nirvana_path}/Data/SupplementaryAnnotation/GRCh38/",
        Nirvana_ref = f"{nirvana_path}/Data/References/Homo_sapiens.GRCh38.Nirvana.dat",
        Nirvana_cache = f"{nirvana_path}/Data/Cache/GRCh38/Both",
        Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll",
        file_name = "05_Annotation/Nirvana/Annotation"

    shell:
        """
        dotnet {params.Nirvana_cmd} \
            -i {input} \
            -o {params.file_name} \
            -r {params.Nirvana_ref} \
            --sd {params.Nirvana_supplementray} \
            -c {params.Nirvana_cache}
        """

rule Annovar:
    input:
        "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output:
        directory("05_Annotation/ANNOVAR")
    threads: 4
    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        annovar_dir = "~/annovar",
        protocol = "refGene,avsnp150,clinvar_20221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,EUR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08",
        operation = "g,f,f,f,f,f,f,f,f,f",
        output = "05_Annotation/ANNOVAR/annotations"

    shell:
        """
        mkdir -p {output}
        perl {params.annovar_dir}/table_annovar.pl {input} {params.annovar_dir}/humandb/ \
            -buildver hg38 \
            -out {params.output} -remove \
            -protocol {params.protocol} \
            -operation {params.operation} \
            -nastring . \
            -vcfinput \
            --thread {threads}
        """

rule Split_variants_idnel:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output: "04_calling/indels/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads:1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include INDEL \
            -O {output}

        """

rule Split_variants_snp:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output: "04_calling/snvs/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads:1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include SNP \
            -O {output}
        """


## to do :
### also add auto sample table (py script)


# rule vep_annotation:
#     input: "04_calling/variants_genotyped.gvcf.gz"
    
    # conda: "GUAP"

#     output: "05_annotation/variants_genotyped.txt"
#     resources:
#         mem_mb=2048,
#         cores=4,
#         mem_gb=2,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt
#     threads:4
#     shell:
#     """
#     vep -i {input} -o {output} --cache \
#         --species homo_sapiens \
#         --assembly GRCh38
#     """

rule multiqc:
    input:
        get_final_output
    
    conda: "wes"

    output:
        "multiqc/multiqc_report.html"
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        "multiqc . -o multiqc/"
