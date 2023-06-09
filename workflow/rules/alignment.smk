

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
    input: f"{ref_fasta_path}/hg38.fa"
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

    output: directory(f"{ref_bwa_path}/{ref_prefix}")
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

    
    conda: "env/wes_gatk.yml"

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

    
    conda: "env/wes_gatk.yml"

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

    
    conda: "env/wes_gatk.yml"

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

