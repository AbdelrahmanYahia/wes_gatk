
rule VariantEval:
    input: "04_calling/{sample}_raw.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/QC/{sample}.eval.grp"
    threads: 1
    params:
        ref = ref_fasta
    benchmark: "benchamrks/VariantEval/{sample}/rawvcf.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        touch {output}
        gatk VariantEval \
            -R {params.ref} \
            --eval:{wildcards.sample} {input} \
            -O {output}
        """

# rule variant_filteration:
#     input: "04_calling/variants_genotyped.gvcf.gz"
    
#     conda: "../env/wes_gatk.yml"

#     output: "04_calling/variants_genotyped_filtered.gvcf.gz"
#     params:
#         ref = ref_fasta
#     benchmark: "benchamrks/variant_filteration/genotypedvcf.txt"

#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     threads: 2
#     shell:
#         """
#         gatk VariantFiltration \
#             --variant {input} \
#             --filter-expression "QD < 2.0"              --filter-name "QD2" \
#             --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
#             --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
#             --filter-expression "FS > 60.0"             --filter-name "FS60" \
#             --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
#             --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
#             --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
#             --create-output-variant-index true \
#             --output {output}
#         """

rule Split_variants_idnel:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Split_variants_idnel/genotypedvcf.txt"

    output: "04_calling/indels/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include INDEL \
            -O {output}

        """

rule variant_filteration_indels:
    input: "04_calling/indels/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/variant_filteration_indels/indels.txt"

    output: "04_calling/indels/variants_genotyped.filttered.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk VariantFiltration \
            --variant {input} \
            --filter-expression "QD < 2.0"                  --filter-name "QD2" \
            --filter-expression "QUAL < 30.0"               --filter-name "QUAL30" \
            --filter-expression "FS > 200.0"                --filter-name "FS200" \
            --filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
            --create-output-variant-index true \
            --output {output}
        """


rule Split_variants_snp:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Split_variants_snp/genotypedvcf.txt"

    output: "04_calling/snvs/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include SNP \
            -O {output}
        """

rule variant_filteration_snps:
    input: "04_calling/snvs/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/variant_filteration_snps/snps.txt"

    output: "04_calling/snvs/variants_genotyped.filttered.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk VariantFiltration \
            --variant {input} \
            --filter-expression "QD < 2.0"              --filter-name "QD2" \
            --filter-expression "QUAL < 30.0"           --filter-name "QUAL30" \
            --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
            --filter-expression "FS > 60.0"             --filter-name "FS60" \
            --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            --create-output-variant-index true \
            --output {output}
        """