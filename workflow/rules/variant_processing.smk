
rule VariantEval:
    input: "04_calling/{sample}_raw.gvcf.gz"
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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


rule Split_variants_idnel:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "env/wes_gatk.yml"

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
    
    conda: "env/wes_gatk.yml"

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

