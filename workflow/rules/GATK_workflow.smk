
rule FqtoUBam:
    input:
        R1 = "00_trimmomatic/{sample}/{sample}_{unit}_1.trimmed.fastq",
        R2 = "00_trimmomatic/{sample}/{sample}_{unit}_2.trimmed.fastq"
    
    conda: "../env/wes_gatk.yml"

    output:
        ubam = "0_samples/{sample}/{sample}_{unit}.ubam"

    benchmark: "benchamrks/FastqToSam/{sample}/{sample}_{unit}.txt"
    threads: 4

    resources:
        mem_mb=5120,
        cores=4,
        mem_gb=5,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    log: 
        "logs/trimmomatic/{sample}/{sample}_{unit}.txt"
    params:
        extra_args = ""
    shell:
        """
        R1={input.R1}
        SM={wildcards.sample}
        LB="{wildcards.sample}_{wildcards.unit}"
        RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)

        picard FastqToSam \
            -F1 {input.R1} -F2 {input.R2} \
            -O {output.ubam} \
            -SM $SM \
            -PL  "illumina" \
            -RG $RGID \
            -LB $LB \
            {params.extra_args}
        """

rule merge_bams:
    input:
        ubam = "0_samples/{sample}/{sample}_{unit}.ubam",
        alignedBam = "02_alignment/{sample}/{sample}_{unit}.bam"    
    conda: "../env/wes_gatk.yml"

    output:
        temp("02_alignment/{sample}/{sample}_{unit}_mergedUnmapped.bam")

    threads: 4
    params:
        fa = ref_fasta

    resources:
        mem_mb=32768,
        cores=4,
        mem_gb=32,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            FastqToSam \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ALIGNED_BAM {input.alignedBam} \
            --UNMAPPED_BAM {input.ubam} \
            --OUTPUT {output} \
            --REFERENCE_SEQUENCE {params.fa} \
            --PAIRED_RUN true \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --CLIP_ADAPTERS false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --ADD_MATE_CIGAR true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --PROGRAM_RECORD_ID "bwamem" \
            --PROGRAM_GROUP_NAME "bwamem" \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true
        """


rule SortNFix:
    input:
        "02_alignment/{sample}/{sample}_{unit}_mergedUnmapped.bam"  
    conda: "../env/wes_gatk.yml"
    output:
        "02_alignment/{sample}/{sample}_{unit}_mergedUnmapped_sorted.bam"

    threads: 4
    params:
        fa = ref_fasta

    resources:
        mem_mb=32768,
        cores=4,
        mem_gb=32,
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
