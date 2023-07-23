rule PreprocessIntervals:
    input:
        ref = ref_fasta
    output:
        "000_CNV_Intervalsprep/targets.preprocessed.interval_list"
    params:
        bed = bed_file

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        mkdir -p 000_CNV_Intervalsprep
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            PreprocessIntervals \
            -R {input.ref} \
            -L {params.bed} \
            --padding 250 \
            --bin-length 0 \
            -imr OVERLAPPING_ONLY \
            -O {output}
        """

rule CollectReadCounts:
    input:
        intervals = "000_CNV_Intervalsprep/targets.preprocessed.interval_list",
        sample = "03_bamPrep/{sample}.pqsr.bam"
    output:
        "06_cnvPrep/{sample}.pqsr.hdf5"
    params:
        ref = ref_fasta,
        bed = bed_file,

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CollectReadCounts \
            -R {params.ref} \
            -L {input.intervals} \
            -imr OVERLAPPING_ONLY \
            -I {input.sample} \
            -O {output}
        """
    
rule AnnotateIntervals:
    input:
        intervals = "000_CNV_Intervalsprep/targets.preprocessed.interval_list"
    output:
        "000_CNV_Intervalsprep/annotated.ref.tsv"
    params:
        ref = ref_fasta

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            AnnotateIntervals \
            -L {input.intervals} \
            -R {params.ref} \
            -imr OVERLAPPING_ONLY \
            -O {output}
        """

rule FilterIntervals:
    input:
        annotatedIntervals = "000_CNV_Intervalsprep/annotated.ref.tsv",
        preprocessedIntervals = "000_CNV_Intervalsprep/targets.preprocessed.interval_list",
        all_samples = expand("06_cnvPrep/{sample}.pqsr.hdf5", sample=set(samples_IDs))

    output:
        "000_CNV_Intervalsprep/ref.cohort.gc.filtered.interval_list"

    params:
        ref = ref_fasta,
        samples = lambda wildcards: [f" -I 06_cnvPrep/{sample}.pqsr.hdf5" for sample in set(samples_IDs)]

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        echo {params.samples} 
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            FilterIntervals \
            -L {input.preprocessedIntervals}\
            --annotated-intervals {input.annotatedIntervals} \
            {params.samples} \
            -imr OVERLAPPING_ONLY \
            --minimum-gc-content 0.1 \
            --maximum-gc-content 0.9 \
            --minimum-mappability 0.9 \
            --maximum-mappability 1.0 \
            --minimum-segmental-duplication-content 0.0 \
            --maximum-segmental-duplication-content 0.5 \
            --low-count-filter-count-threshold 5 \
            --low-count-filter-percentage-of-samples 90.0 \
            --extreme-count-filter-minimum-percentile 1.0 \
            --extreme-count-filter-maximum-percentile 99.0 \
            --extreme-count-filter-percentage-of-samples 90.0 \
            -O {output}
        """

rule DetermineGermlineContigPloidy:
    input:
        filteredIntervals = "000_CNV_Intervalsprep/ref.cohort.gc.filtered.interval_list",
        all_samples = expand("06_cnvPrep/{sample}.pqsr.hdf5", sample=samples_IDs)
        # TODO: double check which samples to include?

    output:
        dir = directory("06_cnvPrep/ploidy-priors"),
        for_automation = "06_cnvPrep/ploidy-priors/ploidy-calls"

    params:
        samples = lambda wildcards: [f" -I 06_cnvPrep/{sample}.pqsr.hdf5" for sample in samples_IDs],
        contigPloidyPriors = config["contig_ploidy_priors"],
        prefix = "06_cnvPrep/ploidy-priors/ploidy" 

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # eval "$(conda shell.bash hook)"
        # set +u; conda activate gatk ; set -u

        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            DetermineGermlineContigPloidy  \
            -L {input.filteredIntervals}\
            {params.samples}  \
            -imr OVERLAPPING_ONLY \
            --contig-ploidy-priors {params.contigPloidyPriors} \
            -O {output.dir} \
            --output-prefix ploidy \
            --verbosity DEBUG
        """

rule GermlineCNVCaller:
    input:
        sample = expand("06_cnvPrep/{sample}.pqsr.hdf5", sample=samples_IDs),
        annotated = "000_CNV_Intervalsprep/annotated.ref.tsv",
        ploidy = '06_cnvPrep/ploidy-priors/ploidy-calls',
        interval = "000_CNV_Intervalsprep/ref.cohort.gc.filtered.interval_list"

    output:
        modelf = "07_cnv/cohort-calls/Allsamples-model",
        callsf = "07_cnv/cohort-calls/Allsamples-calls"

    params:
        outdir = '07_cnv/cohort-calls',
        outpref = 'Allsamples',
        samples = lambda wildcards: [f" -I 06_cnvPrep/{sample}.pqsr.hdf5" for sample in samples_IDs]

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # eval "$(conda shell.bash hook)"
        # set +u; conda activate gatk; set -u

        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GermlineCNVCaller  \
            --run-mode COHORT \
            -L {input.interval} \
            -I {params.samples} \
            --contig-ploidy-calls {input.ploidy}/dogs-calls \
            --annotated-intervals {input.annotated} \
            --output-prefix {params.outpref} \
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {params.outdir}
        """

