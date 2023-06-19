
rule gather_reports:
    input:
        lambda wildcards: expand(
            "03_bamPrep/{sample}/{sample}_{unit}.report",
            unit=units.loc[wildcards.sample, "unit"].tolist(),
            sample=wildcards.sample
        )
    
    conda: "../env/wes_gatk.yml"

    output:
        "03_bamPrep/PQSR.report"
    params:
        reports = lambda wildcards: [f" -I 03_bamPrep/{wildcards.sample}/{wildcards.sample}_{b}.report" for b in units.loc[wildcards.sample, "unit"].tolist()],
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
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GatherBQSRReports \
            {params.reprots} \
            -O {output}
        """
