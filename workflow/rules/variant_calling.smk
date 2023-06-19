
rule HaplotypeCaller:
    input: "03_bamPrep/merged_bams/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        extra_args = config["caller_extra_args"],
        padding = config["padding"]

    resources:
        mem_mb=int(config["calling_mem"])*1024,
        cores=config["calling_threads"],
        mem_gb=config["calling_mem"],
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    threads: config["calling_threads"]

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller -R {params.ref} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation {params.extra_args} \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output}
        """

# same as GenomicsDBImport
rule combine_gvcf:
    input:
        gvcfs=get_gvcf
    
    conda: "../env/wes_gatk.yml"

    output:
        "04_calling/variants.gvcf.gz"
    params:
        gvcfs = lambda wildcards, input: [f" --variant {v}" for v in input["gvcfs"]],
        ref = ref_fasta

    threads: config["general_low_threads"]
    resources:
        mem_mb=int(config["general_low_mem"])* 1024,
        cores=config["general_low_threads"],
        mem_gb=int(config["general_low_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        echo "{params.gvcfs}"
        gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
            CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output}
        """

rule genotype_gvcfs:
    input:
        "04_calling/variants.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        "04_calling/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    threads: config["general_high_threads"]
    resources:
        mem_mb=int(config["general_high_mem"])* 1024,
        cores=config["general_high_threads"],
        mem_gb=int(config["general_high_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GenotypeGVCFs \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            -R {params.ref} -V {input} -O {output}
        """
