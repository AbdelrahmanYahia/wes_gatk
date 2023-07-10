
rule MergeSamFiles:
    input:
        get_merge_input
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Mergesamfiles/{sample}.txt"

    output:
        "03_bamPrep/{sample}.bam"
    params:
        bams = lambda wildcards: [f" -I 02_alignment/{wildcards.sample}/{wildcards.sample}-{b}_mergedUnmapped.bam" for b in units.loc[wildcards.sample, "unit"].tolist()],
        ref = ref_fasta

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        picard MergeSamFiles  \
            {params.bams} \
            -OUTPUT {output} \
            --SORT_ORDER "unsorted" 
        """


rule MarkDuplicates:
    input:
        "03_bamPrep/{sample}.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        bam = "03_bamPrep/{sample}.dedub.bam",
        matrix = "03_bamPrep/{sample}.dedub.matrix"

    threads: 4
    benchmark: "benchamrks/mrkDuplicates/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    log: 
        "logs/mrk-dub/{sample}.dedub.log"

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true > {log}
        """


rule SortNFix:
    input:
        "03_bamPrep/{sample}.dedub.bam"
        # "03_bamPrep/{sample}_mergedUnmapped.bam"  
    conda: "../env/wes_gatk.yml"
    output:
        "03_bamPrep/{sample}.dedub.sorted.bam"

    threads: 4
    params:
        fa = ref_fasta

    benchmark: "benchamrks/SortNFix/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
    gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
      SortSam \
      --INPUT {input} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false | gatk \
      --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE {params.fa}
        """


rule BaseRecalibrator:
    input: "03_bamPrep/{sample}.dedub.sorted.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}.report"
    params: 
        known_sites = known_variants_snps,
        known_sites2 = known_variants_indels,
        known_sites3 = known_variants_indels2,
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]

    threads: 4
    benchmark: "benchamrks/BaseRecalibrator/{sample}.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            BaseRecalibrator \
            -R {params.ref} \
            -I {input} \
            --use-original-qualities \
            --known-sites {params.known_sites} \
            --known-sites {params.known_sites2} \
            --known-sites {params.known_sites3} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -O {output} 
        """

rule applyBaseRecalibrator:
    input: 
        bam = "03_bamPrep/{sample}.dedub.sorted.bam",
        report = "03_bamPrep/{sample}.report"
    benchmark: "benchamrks/ApplyBaseRecab/{sample}.txt"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}.pqsr.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    params: 
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ApplyBQSR  -R {params.ref} \
            -I {input.bam} \
            -bqsr {input.report} -O {output} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 \
            --use-original-qualities
        """


rule bqsr_calibrated_report:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/recaliprated_report/{sample}.txt"

    output: "03_bamPrep/{sample}_pqsr.report"
    params: 
        known_sites = known_variants_snps,
        known_sites2 = known_variants_indels,
        known_sites3 = known_variants_indels2,
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" BaseRecalibrator -R {params.ref} \
            -I {input} \
            --known-sites {params.known_sites} \
            --known-sites {params.known_sites2} \
            --known-sites {params.known_sites3} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -O {output} 
        """

rule AnalyzeCovariates:
    input: 
        raw = "03_bamPrep/{sample}.report", 
        bqsr = "03_bamPrep/{sample}_pqsr.report"

    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/analyzeCovariates/{sample}.txt"

    output: "03_bamPrep/QC/{sample}.pdf"

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" AnalyzeCovariates -before {input.raw} \
            -after {input.bqsr} -plots {output}
        """


# ## TODO: check difference between this and above
