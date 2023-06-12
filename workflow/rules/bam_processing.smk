
rule mrk_duplicates:
    input:
        "02_alignment/{sample}_{unit}.sorted.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        bam = "02_alignment/{sample}_{unit}.dedub.bam",
        matrix = "02_alignment/{sample}_{unit}.dedub.matrix"

    resources:
        mem_mb=4096,
        cores=4,
        mem_gb=4,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    log: 
        "logs/{sample}_{unit}.dedub.log"

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} > {log}
        """

rule BaseRecalibrator:
    input: "02_alignment/{sample}_{unit}.dedub.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}_{unit}.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta,

    threads: 4
    benchmark: "benchamrks/{sample}_{unit}_GATK_pqsr.txt"
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
        bam = "02_alignment/{sample}_{unit}.dedub.bam",
        report = "03_bamPrep/{sample}_{unit}.report"
    benchmark: "benchamrks/{sample}_{unit}_GATK_apply_BQSR.txt"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}_{unit}.pqsr.bam"
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
    input: "03_bamPrep/{sample}_{unit}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}_{unit}_pqsr.report"
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
        raw = "03_bamPrep/{sample}_{unit}.report", 
        bqsr = "03_bamPrep/{sample}_{unit}_pqsr.report"

    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/QC/{sample}_{unit}.pdf"

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
        "03_bamPrep/{sample}_{unit}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("03_bamPrep/QC/{sample}_{unit}_Qualimap")

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

rule MergeSamFiles:
    input:
        get_merge_input
    
    conda: "../env/wes_gatk.yml"

    output:
        "03_bamPrep/merged_bams/{sample}.pqsr.bam"
    params:
        bams = lambda wildcards: [f" -I 03_bamPrep/{wildcards.sample}_{b}.pqsr.bam" for b in units.loc[wildcards.sample, "unit"].tolist()],
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
        picard MergeSamFiles  \
            {params.bams} \
            -OUTPUT {output} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true
        """
