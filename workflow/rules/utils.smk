
def get_gvcf(wildcards):
    return [f"04_calling/{sample}_raw.gvcf.gz" \
        for sample in samples_IDs]

# def get_bams(wildcards):
#     return [f"03_bamPrep/{sample}.pqsr.bam" \
#         for sample in samples_IDs]

def get_final_output(wildcards):
    final_output = []

    for i in samples_IDs:
        # final_output.extend(expand(
        #     "01_QC/{sample}/{sample}_{unit}_{R}.trimmed_fastqc.zip",
        #     sample=i,
        #     unit=units.loc[i, "unit"].tolist(),
        #     R=[1,2]
        # ))
        final_output.extend(expand(
            "02_alignment/{sample}/QC/{sample}_{unit}_mergedUnmapped.cov",
                sample = i,
                unit=units.loc[i, "unit"].tolist()
        ))
        final_output.extend(expand(
            "03_bamPrep/{sample}/QC/{sample}_{unit}_Qualimap",
                sample = i,
                unit=units.loc[i, "unit"].tolist()
        ))

        final_output.extend(expand(
            "03_bamPrep/{sample}/QC/{sample}_{unit}.pdf",
                sample = i,
                unit=units.loc[i, "unit"].tolist()
        ))
    
    final_output.extend(expand(
        "04_calling/QC/{sample}.eval.grp",
            sample = samples_IDs
    ))
    final_output.extend(expand(
            "04_calling/QC/{type}/bcftools.stats",
            type = ["snvs", "indels"]
    ))
    final_output.extend(expand(
            "05_Annotation/ANNOVAR/{type}",
            type = ["snvs", "indels"]
    ))
    final_output.extend(expand(
            "05_Annotation/Nirvana/{type}/Annotation.json.gz",
            type = ["snvs", "indels"]
    ))


    return final_output


def get_merge_input(wildcards):

    return expand(
            "03_bamPrep/{sample}/{sample}_{unit}.pqsr.bam",
            unit=units.loc[wildcards.sample, "unit"].tolist(),
            sample=wildcards.sample
        )

