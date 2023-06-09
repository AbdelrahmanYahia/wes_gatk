## to do :
### also add auto sample table (py script)


# rule vep_annotation:
#     input: "04_calling/variants_genotyped.gvcf.gz"
    
    # conda: "GUAP"

#     output: "05_annotation/variants_genotyped.txt"
#     resources:
#         mem_mb=2048,
#         cores=4,
#         mem_gb=2,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt
#     threads:4
#     shell:
#     """
#     vep -i {input} -o {output} --cache \
#         --species homo_sapiens \
#         --assembly GRCh38
#     """

