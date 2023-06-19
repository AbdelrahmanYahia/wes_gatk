## TODO: implememt this workflow
# rule merge_bams:
#     input:
#         ubam = "0_samples/{sample}/{sample}-{unit}.ubam",
#         alignedBam = "02_alignment/{sample}/{sample}-{unit}.bam"    
#     conda: "../env/wes_gatk.yml"

#     output:
#         temp("02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam")

#     threads: 4
#     params:
#         fa = ref_fasta

#     resources:
#         mem_mb=4096,
#         cores=4,
#         mem_gb=4,
#         nodes = 1,
#         time = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             FastqToSam \
#             --VALIDATION_STRINGENCY SILENT \
#             --EXPECTED_ORIENTATIONS FR \
#             --ATTRIBUTES_TO_RETAIN X0 \
#             --ALIGNED_BAM {input.alignedBam} \
#             --UNMAPPED_BAM {input.ubam} \
#             --OUTPUT {output} \
#             --REFERENCE_SEQUENCE {params.fa} \
#             --PAIRED_RUN true \
#             --SORT_ORDER "unsorted" \
#             --IS_BISULFITE_SEQUENCE false \
#             --ALIGNED_READS_ONLY false \
#             --CLIP_ADAPTERS false \
#             --MAX_RECORDS_IN_RAM 2000000 \
#             --ADD_MATE_CIGAR true \
#             --MAX_INSERTIONS_OR_DELETIONS -1 \
#             --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
#             --PROGRAM_RECORD_ID "bwamem" \
#             --PROGRAM_GROUP_NAME "bwamem" \
#             --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
#             --ALIGNER_PROPER_PAIR_FLAGS true \
#             --UNMAP_CONTAMINANT_READS true
#         """

rule SortNFix:
    input:
        "02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam"  
    conda: "../env/wes_gatk.yml"
    output:
        "02_alignment/{sample}/{sample}-{unit}_mergedUnmapped_sorted.bam"

    threads: 4
    params:
        fa = ref_fasta

    resources:
        mem_mb=4096,
        cores=4,
        mem_gb=4,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
    gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
      SortSam \
      --INPUT {input} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false | gatk \
      --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE {params.fa}
        """


rule mrk_duplicates:
    input:
        "02_alignment/{sample}/{sample}-{unit}_mergedUnmapped_sorted.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        bam = "02_alignment/{sample}/{sample}-{unit}.dedub.bam",
        matrix = "02_alignment/{sample}/{sample}-{unit}.dedub.matrix"
    params: 
        extra_args = """--VALIDATION_STRINGENCY SILENT \
                        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                        --ASSUME_SORT_ORDER "queryname" \
                        --CREATE_MD5_FILE true""",

    resources:
        mem_mb=4096,
        cores=4,
        mem_gb=4,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    log: 
        "logs/mrk-dub/{sample}/{sample}-{unit}.dedub.log"

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} \
            {params.extra_args} > {log}
        """

rule BaseRecalibrator:
    input: "02_alignment/{sample}/{sample}-{unit}.dedub.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}/{sample}-{unit}.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta,

    threads: 4
    benchmark: "benchamrks/{sample}/{sample}-{unit}_GATK_pqsr.txt"
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            BaseRecalibrator \
            -R {params.ref} \
            -I {input} \
            --use-original-qualities \
            --known-sites {params.known_sites} \
            -O {output} 
        """

rule applyBaseRecalibrator:
    input: 
        bam = "02_alignment/{sample}/{sample}-{unit}.dedub.bam",
        report = "03_bamPrep/{sample}/{sample}-{unit}.report"
    benchmark: "benchamrks/{sample}/{sample}-{unit}_GATK_apply_BQSR.txt"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}/{sample}-{unit}.pqsr.bam"
    threads: 1
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ApplyBQSR  -R {params.ref} \
            -I {input.bam} --emit-original-quals \
            -bqsr {input.report} -O {output} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 \
            --use-original-qualities
        """


rule bqsr_calibrated_report:
    input: "03_bamPrep/{sample}/{sample}-{unit}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}/{sample}-{unit}_pqsr.report"
    params: 
        known_sites = known_variants,
        ref = ref_fasta
    threads: 1
    resources:
        mem_mb=2048,
        cores=4,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" BaseRecalibrator -R {params.ref} \
            -I {input} --known-sites {params.known_sites} \
            -O {output} 
        """

rule AnalyzeCovariates:
    input: 
        raw = "03_bamPrep/{sample}/{sample}-{unit}.report", 
        bqsr = "03_bamPrep/{sample}/{sample}-{unit}_pqsr.report"

    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}/QC/{sample}-{unit}.pdf"

    threads: 2
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" AnalyzeCovariates -before {input.raw} \
            -after {input.bqsr} -plots {output}
        """

rule qualimap:
    input:
        "03_bamPrep/{sample}/{sample}-{unit}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("03_bamPrep/{sample}/QC/{sample}-{unit}_Qualimap")

    threads: 2
    resources:
        mem_mb=2048,
        cores=2,
        mem_gb=2,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt  

    shell:
        """
        qualimap \
            bamqc \
            -bam {input} \
            --paint-chromosome-limits \
            --genome-gc-distr HUMAN \
            -nt {threads} \
            -skip-duplicated \
            --skip-dup-mode 0 \
            -outdir {output} \
            -outformat HTML
        """


rule MergeSamFiles:
    input:
        get_merge_input
    
    conda: "../env/wes_gatk.yml"

    output:
        "03_bamPrep/merged_bams/{sample}.pqsr.bam"
    params:
        bams = lambda wildcards: [f" -I 03_bamPrep/{wildcards.sample}/{wildcards.sample}-{b}.pqsr.bam" for b in units.loc[wildcards.sample, "unit"].tolist()],
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
        picard MergeSamFiles  \
            {params.bams} \
            -OUTPUT {output} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true
        """

## TODO: check difference between this and above

rule GatherBamFiles:
    input:
        get_merge_input
    
    conda: "../env/wes_gatk.yml"

    output:
        "03_bamPrep/merged_bams/{sample}.gatherd.pqsr.bam"
    params:
        bams = lambda wildcards: [f" -I 03_bamPrep/{wildcards.sample}/{wildcards.sample}-{b}.pqsr.bam" for b in units.loc[wildcards.sample, "unit"].tolist()],
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
            GatherBamFiles \
            {params.bams} \
            -OUTPUT {output} \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true
        """

