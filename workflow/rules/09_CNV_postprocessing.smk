
def sampleindex(sample):
    samples = list(samples_IDs)
    index = samples.index(sample)
    return index

rule PostprocessGermlineCNVCalls:
    input:
        model = "07_cnv/cohort-calls/Allsamples-model",
        calls = "07_cnv/cohort-calls/Allsamples-calls",
        seq_dict = str(config["reference_fasta"]).replace(".fa", ".dict")

    output:
        intervals = "08_cnv_postprocessing/{sample}.intervals_cohort.vcf.gz",
        segments = "08_cnv_postprocessing/{sample}.segments_cohort.vcf.gz"
    params:
        index = lambda wildcards: sampleindex(wildcards.sample)

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            PostprocessGermlineCNVCalls \
            --model-shard-path {input.model} \
            --calls-shard-path {input.call} \
            --allosomal-contig chrX --allosomal-contig chrY \
            --contig-ploidy-calls ploidy-calls \
            --sample-index {params.index} \
            --output-genotyped-intervals  {output.intervals} \
            --output-genotyped-segments  {output.segments} \
            --sequence-dictionary {input.seq_dict}

        """

## you can visualize the gCNV with IGV


### you can perfrom the following:
# [i] VariantsToTable to subset and columnize annotations

rule VariantsToTable :
    input:
        "08_cnv_postprocessing/{sample}.segments_cohort.vcf.gz",

    output:
        "08_cnv_postprocessing/{sample}/genotyped-segments-case-{sample}-vs-cohort.table.txt",

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            VariantsToTable \
            -V {input} \
            -F CHROM -F POS -F END -GF NP -GF CN \
            -O {output}

        """

# [ii] Unix shell commands to convert to SEG format data
rule VariantsToTable_seg:
    input:
        txt = "08_cnv_postprocessing/{sample}/genotyped-segments-case-{sample}-vs-cohort.table.txt",
        vcf = "08_cnv_postprocessing/{sample}.segments_cohort.vcf.gz"
    output:
        "08_cnv_postprocessing/{sample}/genotyped-segments-case-{sample}-vs-cohort.table.txt.seg"

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        
        sampleName=$(gzcat {input.vcf} | grep -v '##' | head -n1 | cut -f10)
        awk -v sampleName=$sampleName 'BEGIN {FS=OFS="\t"} {print sampleName, $0}' {input.txt} &gt; ${i}.seg; head ${i}.seg

        """
