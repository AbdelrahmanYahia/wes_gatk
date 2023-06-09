
def get_gvcf(wildcards):
    return [f"04_calling/{sample}_raw.gvcf.gz" \
        for sample in samples]

def get_bams(wildcards):
    return [f"03_bamPrep/{sample}.pqsr.bam" \
        for sample in samples]

def get_final_output(wildcards):
    final_output = []
    final_output.extend(expand(
        f"00_QC/{{sample}}_R{{R}}_fastqc.zip",
            sample = samples,
            R = [1,2]
    ))
    final_output.extend(expand(
        "02_alignment/QC/{sample}.cov",
            sample = samples
    ))
    final_output.extend(expand(
        "03_bamPrep/QC/{sample}_Qualimap",
        sample = samples
    ))
    final_output.extend(expand(
        "03_bamPrep/QC/{sample}.pdf",
            sample = samples
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
            sample = samples
    ))



    return final_output
