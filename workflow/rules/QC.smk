
rule trimmomatic:
    input:
        R1 = f"{samples_dir}/{{sample}}_R1.fastq",
        R2 = f"{samples_dir}/{{sample}}_R2.fastq"
    
    conda: "env/wes_gatk.yml"

    output:
        log="logs/trimmomatic/{sample}_{unit}.log",
        summary="logs/trimmomatic/{sample}.summary",
        nf1 = "00_trimmomatic/{sample}_{unit}_1.trimmed.fastq.gz",
        nf2 = "00_trimmomatic/{sample}_{unit}_2.trimmed.fastq.gz",
        nfu1=temp("00_trimmomatic/{sample}_{unit}-U_1.trimmed.fastq.gz"),
        nfu2=temp("00_trimmomatic/{sample}_{unit}-U_1.trimmed.fastq.gz")

    benchmark: "benchamrks/QC/{sample}_{unit}_trim.txt"
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
        "logs/trimmomatic/{sample}_{unit}.txt"
    shell:
        """
        trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
                -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
                SLIDINGWINDOW:{params.size}:{params.quality} MINLEN:{params.minlen} > {log} 2>&1
        """

rule Fastqc:
    input:
        "00_trimmomatic/{sample}_{unit}_{R}.trimmed.fastq.gz"
    log:
        "logs/QC/QC_{sample}_{unit}_{R}.log"
    
    conda: "env/wes_gatk.yml"

    output:
        zip="01_QC/{sample}_{unit}_{R}_fastqc.zip",
        html="01_QC/{sample}_{unit}_{R}_fastqc.html"
        
    benchmark: f"benchamrks/QC/{sample}_{unit}_{R}_fastqc.txt"
    threads: 2
    params:
        path="01_QC"
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        "fastqc {input} --threads {threads} -o {params.path} > {log} 2>&1"


rule multiqc:
    input:
        get_final_output
    
    conda: "env/wes_gatk.yml"

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