
rule HaplotypeCaller:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        extra_args = config["caller_extra_args"],
        padding = config["padding"]

    benchmark: "benchamrks/HaplotypeCaller/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller -R {params.ref} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation {params.extra_args} \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output}
        """

rule ReblockGVCF:
    input: "04_calling/{sample}_raw.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04-1_gvcf-processing/{sample}_reblocked.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]

    benchmark: "benchamrks/Reblock/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    threads: 4

    shell:
        """
        ## IMPORTANT!
        # Note that when uncalled alleles are dropped, 
        # the original GQ may increase. Use --keep-all-alts 
        # if GQ accuracy is a concern.


        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ReblockGVCF \
            -R {params.ref} \
            -G StandardAnnotation \
            -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            --floor-blocks \
            -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -V {input} -O {output}
        """

rule CreateSampleSheet:
    input:
        gvcfs=get_reblocked_gvcf
    
    conda: "../env/wes_gatk.yml"

    output:
        "resources/sample_sheet.tsv"

    params:
        IDs = samples_IDs

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        touch {output}

        for i in {params.IDs} ; do
            echo -e  "${{i}}\t04-1_gvcf-processing/${{i}}_reblocked.gvcf.gz" | cat >> {output}
        done
        """

rule scatter_intervals:
    input:
        ref = ref_fasta
    
    conda: "../env/wes_gatk.yml"

    output:
        dir = directory("resources/Scatterred_intervals"),
        intervals = expand("resources/Scatterred_intervals/00{interval}-scattered.interval_list", interval = [f"{i:02d}" for i in range(0, int(config['scatter_count']))])

    params:
        scatter_count = int(config['scatter_count']),
        use_interval = lambda wildcards: f"-L {bed_file}" if config['use_supplied_interval'] else "",
        use_interval = lambda wildcards: [f" -L {v}" for v in ],
    
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

        # gatk SplitIntervals -L {input.intervals} \
        # -O {output.out_dir} -R {input.ref} -scatter {params} \
        # -mode INTERVAL_SUBDIVISION \
        # --interval-merging-rule OVERLAPPING_ONLY &> {log}
        # ls -l {output.out_dir}/*scattered.interval_list > {output.fof}
        # --subdivision-mode,-mode <IntervalListScatterMode>
        #                       How to divide intervals.  Default value: INTERVAL_SUBDIVISION. INTERVAL_SUBDIVISION
        #                       (Scatter the interval list into similarly sized interval lists (by base count), breaking
        #                       up intervals as needed.)
        #                       BALANCING_WITHOUT_INTERVAL_SUBDIVISION (Scatter the interval list into similarly sized
        #                       interval lists (by base count), but without breaking up intervals.)
        #                       BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW (Scatter the interval list into
        #                       similarly sized interval lists (by base count), but without breaking up intervals. Will
        #                       overflow current interval list so that the remaining lists will not have too many bases to
        #                       deal with.)
        #                       INTERVAL_COUNT (Scatter the interval list into similarly sized interval lists (by interval
        #                       count, not by base count). Resulting interval lists will contain the same number of
        #                       intervals except for the last, which contains the remainder.)
        #                       INTERVAL_COUNT_WITH_DISTRIBUTED_REMAINDER (Scatter the interval list into similarly sized
        #                       interval lists (by interval count, not by base count). Resulting interval lists will
        #                       contain similar number of intervals.)
        # --interval-merging-rule OVERLAPPING_ONLY or ALL [default all]
        
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            SplitIntervals \
            -R {input} {params.use_interval} \
            --scatter-count {params.scatter_count} \
            --interval-merging-rule ALL \
            -mode INTERVAL_SUBDIVISION \
            -O {output.dir} 
        """




"""
sort -k1,1 -k2,2n $bedfile | awk '
{
    # If it is a different chromosome or if it is a non-overlapping interval
    if ($1 != chrom) {
        # If it is not the first interval
        if (chrom) {
            print chrom, start, end;
        }
        chrom = $1;
        start = $2;
        end = $3;
    } else if ($3 > end) {
        # If it is an overlapping interval
        end = $3;
    }
}
END {
    # Print the last interval
    print chrom, start, end;
}
' > output.bed
"""





rule GenomicsDBImport:
    input:
        sample_map = "resources/sample_sheet.tsv",
        intervals = "resources/Scatterred_intervals/00{interval}-scattered.interval_list"

    conda: "../env/wes_gatk.yml"

    output:
        directory("04-2_GenomicsDB/Interval_{interval}")

    params:
        ref = ref_fasta,
        interval_optimzation = lambda wildcards: f"--merge-contigs-into-num-partitions 50" if config['use_supplied_interval'] else "",

    benchmark: "benchamrks/GenomicsDBImport/Interval_{interval}.txt"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (32 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 32  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        reduced = lambda wildcards, attempt: int(attempt * (32 * 0.80))
    shell:
        """

        # TODO: we can implement the update here: by checking and replacing:
        ## --genomicsdb-update-workspace-path path_to_dir instead of:
        ## --genomicsdb-workspace-path path_to_dir
        ## but when updating I don't need to add -L as it figures that out from the samples automatically

        # sources: https://eriqande.github.io/eca-bioinf-handbook/variant-calling.html#genomics-db
        # finally  some notes:
        ## the documentation for GenomicsDBImport is replete with warnings 
        ## about how you should backup your genomics data bases before trying
        ## to update them. Apparently there is a good chance that any little
        ## glitch during the updating process could corrupt the entire data base and render it useless. Joy!

        ## With that as a preamble, it does seem that we should try creating some data bases and then adding samples to them.


        ### use --genomicsdb-shared-posixfs-optimizations  for HPC


        mkdir -p 04-2_GenomicsDB
        
        
        gatk --java-options "-Xmx{resources.reduced}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads} -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
            GenomicsDBImport \
            --sample-name-map {input.sample_map} \
            -L {input.intervals} \
            --genomicsdb-workspace-path {output} \
            --batch-size 50 \
            --reader-threads 4 \
            --merge-input-intervals \
            --consolidate \
            --genomicsdb-shared-posixfs-optimizations true

        """
### A USER ERROR has occurred: Input variant files must be block compressed vcfs when using bypass-feature-reader
## when using --bypass-feature-reader will try and fix it later 

rule genotype_gvcfs:
    input:
        "04-2_GenomicsDB/Interval_{interval}"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

    output:
        "04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz"
    params:
        ref = ref_fasta
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GenotypeGVCFs \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            -R {params.ref} -V gendb://{input} -O {output}
        """

## from: https://github.com/harvardinformatics/snpArcher/blob/main/workflow/rules/bam2vcf_gatk_intervals.smk
# rule sort_gatherVcfs:
#     input:
#         vcfs = expand("04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz", interval = [f"{i:02d}" for i in range(0, int(config['scatter_count']))]),
    
#     output:
#         vcfFinal = "04_calling/variants_genotyped.gvcf.gz",
#         vcfFinalidx = "04_calling/variants_genotyped.gvcf.gz.tbi"

#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/sort_gatherVcfs/all_samples.txt"

#     shell:
#         """
#         bcftools concat -D -a -Ou {input.vcfs} | bcftools sort -Oz -o {output.vcfFinal} - 
#         tabix -p vcf {output.vcfFinal} 
#         """

rule combine_gvcf:
    input:
        gvcfs= expand("04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz", interval = [f"{i:02d}" for i in range(0, int(config['scatter_count']))])
    
    conda: "../env/wes_gatk.yml"

    output:
        "04_calling/variants.gvcf.gz"
    params:
        gvcfs = lambda wildcards, input: [f" --variant {v}" for v in input["gvcfs"]],
        ref = ref_fasta

    benchmark: "benchamrks/combine_gvcf/combinedvcf.txt"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        echo "{params.gvcfs}"
        gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
            CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output}
        """
