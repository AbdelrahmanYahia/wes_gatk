
## TODO: add gatk best practice args:
#           -R ~{ref_fasta} \
#           -I ~{input_bam} \
#           -L ~{interval_list} \
#           -O ~{output_filename} \
#           -contamination ~{default="0" contamination} \
#           -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
#           -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
#           ~{true="-ERC GVCF" false="" make_gvcf} \
#           ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
#           ~{bamout_arg}

## TODO: update 
rule HaplotypeCaller:
    input: "03_bamPrep/merged_bams/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        workflow_args = """"""

    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller -R {params.ref} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \

            -L {params.bed} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output}
        """

## TODO: merging samples lanes and libraries?
## TODO: Genotype each sample individualy?

rule combine_gvcf:
    input:
        gvcfs=get_gvcf
    
    conda: "../env/wes_gatk.yml"

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
        gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
            CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
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
    threads: 4
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GenotypeGVCFs \
            -R {params.ref} -V {input} -O {output}
        """
