
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