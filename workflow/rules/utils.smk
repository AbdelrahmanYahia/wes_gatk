
def get_raw_gvcf(wildcards):
    return [f"04_calling/{sample}_raw.gvcf.gz" \
        for sample in set(samples_IDs)]

def get_reblocked_gvcf(wildcards):
    return [f"04-1_gvcf-processing/{sample}_reblocked.gvcf.gz" \
        for sample in set(samples_IDs)]

# def get_bams(wildcards):
#     return [f"03_bamPrep/{sample}.pqsr.bam" \
#         for sample in samples_IDs]

def get_final_output(wildcards):
    final_output = []
    final_output.extend(expand(
            "04_calling/QC/{sample}.eval.grp",
            sample = samples_IDs
    ))

    final_output.extend(expand(
            "05_Annotation/ANNOVAR/{type}",
            type = ["snvs", "indels"]
    ))

    final_output.extend(expand(
            "05_Annotation/Nirvana/{type}/Annotation.json.gz",
            type = ["snvs", "indels"]
    ))

    final_output.extend(expand(
            "03_bamPrep/QC/{sample}.cov",
            sample = samples_IDs
    ))

    final_output.extend(expand(
            "03_bamPrep/QC/{sample}_Qualimap",
            sample = samples_IDs
    ))

    final_output.extend(expand(
            "03_bamPrep/QC/{sample}.pdf",
            sample = samples_IDs
    ))

    if call_cnv:
        # final_output.extend(["07_cnv/cohort-calls/Allsamples-calls"])
        final_output.extend(expand(
            "08_cnv_postprocessing/{sample}.intervals_cohort.vcf.gz",
            sample = samples_IDs
    ))

    # final_output.extend(expand(
    #         "04_calling/QC/{type}/bcftools.stats",
    #         type = ["snvs", "indels"]
    # ))

    # for i in samples_IDs:
        # final_output.extend(expand(
        #     "01_QC/{sample}/{sample}-{unit}_{R}.trimmed_fastqc.zip",
        #     sample=i,
        #     unit=units.loc[i, "unit"].tolist(),
        #     R=[1,2]
        # ))
        

    return final_output


def get_merge_input(wildcards):
    return expand(
            "02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam",
            unit=units.loc[wildcards.sample, "unit"].tolist(),
            sample=wildcards.sample
        )

