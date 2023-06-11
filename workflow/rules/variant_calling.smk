
rule HaplotypeCaller:
    input: "03_bamPrep/merged_bams/{sample}.pqsr.bam"
    
    conda: "env/wes_gatk.yml"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file

    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" HaplotypeCaller -R {params.ref} \
            -L {params.bed} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output}
        """

rule combine_gvcf:
    input:
        gvcfs=get_gvcf
    
    conda: "env/wes_gatk.yml"

    output:
        "04_calling/variants.gvcf.gz"
    params:
        gvcfs = lambda wildcards, input: [f" --variant {v}" for v in input["gvcfs"]],
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
        gatk CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
            -O {output}
        """

rule genotype_gvcfs:
    input:
        "04_calling/variants.gvcf.gz"
    
    conda: "env/wes_gatk.yml"

    output:
        "04_calling/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    threads: 4
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" GenotypeGVCFs \
            -R {params.ref} -V {input} -O {output}
        """
