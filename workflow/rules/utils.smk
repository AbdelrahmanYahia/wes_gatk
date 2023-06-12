
def get_gvcf(wildcards):
    return [f"04_calling/{sample}_raw.gvcf.gz" \
        for sample in samples_IDs]

# def get_bams(wildcards):
#     return [f"03_bamPrep/{sample}.pqsr.bam" \
#         for sample in samples_IDs]

def get_final_output(wildcards):
    final_output = []

    for i in samples_IDs:
        final_output.extend(expand(
            "01_QC/{sample}/{sample}_{unit}_{R}.trimmed_fastqc.zip",
            sample=i,
            unit=units.loc[i, "unit"].tolist(),
            R=[1,2]
        ))
        final_output.extend(expand(
            "02_alignment/{sample}/QC/{sample}_{unit}.cov",
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
    
    final_output.extend(
        [
            "04_calling/QC/bcftools.stats",
            "04_calling/QC/bcftools_plots",
            "04_calling/QC/bcftools_csq.vcf",
            "04_calling/snvs/variants_genotyped.gvcf.gz",
            "04_calling/indels/variants_genotyped.gvcf.gz",
            "05_Annotation/Nirvana/Annotation.json.gz",
            "05_Annotation/ANNOVAR"
        ]
    )
    final_output.extend(expand(
        "04_calling/QC/{sample}.eval.grp",
            sample = samples_IDs
    ))

    return final_output




def get_merge_input(wildcards):

    return expand(
            "03_bamPrep/{sample}/{sample}_{unit}.pqsr.bam",
            unit=units.loc[wildcards.sample, "unit"].tolist(),
            sample=wildcards.sample
        )

