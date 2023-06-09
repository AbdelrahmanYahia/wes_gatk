
rule mrk_duplicates:
    input:
        "02_alignment/{sample}.sorted.bam"
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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

    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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

