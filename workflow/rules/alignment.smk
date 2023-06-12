
rule index_ref:
    input: f"{ref_fasta}"
    
    conda: "env/wes_gatk.yml"

    output: f"{ref_fasta}.fai"

    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: "samtools faidx {input}"


rule refrence_dict:
    input: f"{ref_fasta}"
    
    conda: "env/wes_gatk.yml"

    output: f"{ref_fasta}".replace(".fa", ".dict")
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: "picard CreateSequenceDictionary -R {input}"

rule bwa_index:
    input: f"{ref_fasta}"
    
    conda: "env/wes_gatk.yml"

    output: directory(f"{ref_bwa_path}/{ref_prefix}")
    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell: 
        "bwa index -p {output} {input}"

## lambda <arguments> : <statement1> if <condition> else <statement2>
rule bwa_align:
    input:
        R1 = "00_trimmomatic/{sample}_{unit}_1.trimmed.fastq.gz",
        R2 = "00_trimmomatic/{sample}_{unit}_2.trimmed.fastq.gz"
    
    conda: "env/wes_gatk.yml"

    output:
        temp("02_alignment/{sample}_{unit}.sam")

    threads: 4
    params:
        fa = ref_fasta,
        # library_index = lambda wildcards: units.loc[:, 'library_index'][units['unit'] == f"{wildcards.unit}"].tolist()[0],
        # lane = lambda wildcards: units.loc[:, 'lane'][units['unit'] == f"{wildcards.unit}"].tolist()[0],
        index = lambda wildcards: ref_bwa if not config["index_fasta"] else rules.bwa_index.output

    log: 
        bwa = "logs/{sample}_{unit}_bwa.log",

    benchmark: "benchamrks/{sample}_{unit}_bwa.txt"
    resources:
        mem_mb=32768,
        cores=4,
        mem_gb=32,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        R1={input.R1}
        SM={wildcards.sample}
        PL="Illumina"
        LB={wildcards.unit}
        name=$(basename $R1 | cut -d'_' -f1)
        RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
        PU=$RGID.$LB 
        bwa mem -t {threads} -M \
            -R "@RG\\tID:$RGID\\tSM:$SM\\tPL:$PL\\tLB:$LB\\tPU:$PU" {params.index} {input.R1} {input.R2} > {output} 2> {log.bwa}
        """


rule sort_and_convert_sam:
    input:
        "02_alignment/{sample}_{unit}.sam"

    
    conda: "env/wes_gatk.yml"

    output:
        "02_alignment/{sample}_{unit}.sorted.bam"
    
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
        "02_alignment/{sample}_{unit}.sorted.bam"

    conda: "env/wes_gatk.yml"

    output:
        cov = "02_alignment/QC/{sample}_{unit}.cov",
        stats = "02_alignment/QC/{sample}_{unit}.stats"
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

