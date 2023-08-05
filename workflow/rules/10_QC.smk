
rule multiqc:
    input:
        get_final_output
    
    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/Multiqc/report.txt"

    output:
        "multiqc/multiqc_report.html"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        "multiqc . -o multiqc/"

# rule trimmomatic:
#     input:
#         R1 = f"{samples_dir}/{{sample}}-{{unit}}1.{EXT}",
#         R2 = f"{samples_dir}/{{sample}}-{{unit}}2.{EXT}"
    
#     conda: "../env/wes_gatk.yml"

#     output:
#         log="logs/trimmomatic/{sample}/{sample}-{unit}.log",
#         summary="logs/trimmomatic/{sample}/{sample}-{unit}.summary",
#         nf1 = "00_trimmomatic/{sample}/{sample}-{unit}1.trimmed.fastq.gz",
#         nf2 = "00_trimmomatic/{sample}/{sample}-{unit}2.trimmed.fastq.gz",
#         nfu1=temp("00_trimmomatic/{sample}/{sample}-{unit}-U_1.trimmed.fastq.gz"),
#         nfu2=temp("00_trimmomatic/{sample}/{sample}-{unit}-U_2.trimmed.fastq.gz")

#     benchmark: "benchamrks/QC/{sample}/{sample}-{unit}_trim.txt"
#     threads: 4
#     params:
#         size = 4,
#         quality = 10,
#         extra = '',
#         minlen = 50

#     resources:
#         mem_mb=5120,
#         cores=4,
#         mem_gb=5,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt
#     log: 
#         "logs/trimmomatic/{sample}/{sample}-{unit}.txt"
#     shell:
#         """
#         trimmomatic PE -threads {threads} -phred33 -trimlog {output.log} \
#                 -summary {output.summary} {input.R1} {input.R2} {output.nf1} {output.nfu1} {output.nf2} {output.nfu2} \
#                 SLIDINGWINDOW:{params.size}:{params.quality} MINLEN:{params.minlen} > {log} 2>&1
#         """

# rule gunzip_trimmomatic:
#     input:
#         nf1 = "00_trimmomatic/{sample}/{sample}-{unit}1.trimmed.fastq.gz",
#         nf2 = "00_trimmomatic/{sample}/{sample}-{unit}2.trimmed.fastq.gz",
    
#     conda: "../env/wes_gatk.yml"

#     output:
#         nf1 = "00_trimmomatic/{sample}/{sample}-{unit}1.trimmed.fastq",
#         nf2 = "00_trimmomatic/{sample}/{sample}-{unit}2.trimmed.fastq"

#     threads: 1

#     resources:
#         mem_mb=2048,
#         cores=1,
#         mem_gb=2,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gunzip {input.nf1} 
#         gunzip {input.nf2} 
#         """

# rule FqtoUBam:
#     input:
#         R1 = "00_trimmomatic/{sample}/{sample}-{unit}1.trimmed.fastq",
#         R2 = "00_trimmomatic/{sample}/{sample}-{unit}2.trimmed.fastq"
    
#     conda: "../env/wes_gatk.yml"

#     output:
#         ubam = "0_samples/{sample}/{sample}-{unit}-trimmed.ubam"

#     benchmark: "benchamrks/FastqToSam/{sample}/{sample}-{unit}.txt"
#     threads: 4

#     resources:
#         mem_mb=5120,
#         cores=4,
#         mem_gb=5,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt
#     log: 
#         "logs/trimmomatic/{sample}/{sample}-{unit}.txt"
#     params:
#         extra_args = ""
#     shell:
#         """
#         R1={input.R1}
#         SM={wildcards.sample}
#         LB="{wildcards.sample}_{wildcards.unit}"
#         RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)

#         picard FastqToSam \
#             -F1 {input.R1} -F2 {input.R2} \
#             -O {output.ubam} \
#             -SM $SM \
#             -PL  "illumina" \
#             -RG $RGID \
#             -LB $LB \
#             {params.extra_args}
#         """

# rule Fastqc:
#     input:
#         f"00_trimmomatic/{{sample}}/{{sample}}-{{unit}}{R}.fq.gz"
#     log:
#         "logs/QC/QC_{sample}-{unit}_{R}.log"
    
#     conda: "../env/wes_gatk.yml"

#     output:
#         zip="01_QC/{sample}/{sample}-{unit}{R}.trimmed_fastqc.zip",
#         html="01_QC/{sample}/{sample}-{unit}{R}.trimmed_fastqc.html"
        
#     benchmark: "benchamrks/QC/{sample}/{sample}-{unit}{R}_fastqc.txt"
#     threads: 2
#     params:
#         path=lambda wildcards: f"01_QC/{wildcards.sample}"
#     resources:
#         mem_mb=2048,
#         cores=2,
#         mem_gb=2,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt
#     shell:
#         """
#         mkdir -p {params.path}
#         fastqc {input} --threads {threads} -o {params.path} > {log} 2>&1
#         """
