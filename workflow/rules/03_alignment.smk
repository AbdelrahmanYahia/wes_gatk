
rule ubam_align:
    input:
        bam="0_samples/{sample}/{sample}-{unit}.adab.ubam"

    output:
        bam="02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam"

    conda: "../env/wes_gatk.yml"
    threads: 4
    params:
        fa = ref_fasta,
        index = ref_bwa,
        bwa_args = config["aligner_extra_args"]
        ## TODO: fix this
        # index = lambda wildcards: ref_bwa if not config["index_fasta"] else rules.bwa_index.output

        # library_index = lambda wildcards: units.loc[:, 'library_index'][units['unit'] == f"{wildcards.unit}"].tolist()[0],
        # lane = lambda wildcards: units.loc[:, 'lane'][units['unit'] == f"{wildcards.unit}"].tolist()[0]

    # log: 
    #     bwa = "logs/bwa/{sample}/{sample}-{unit}_gatk-bwa.log",

    benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
    resources:
        mem_mb = 32* 1024,
        # cores=config["align_threads"],
        mem_gb = 32,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            SamToFastq \
            -I {input.bam} \
            --FASTQ /dev/stdout \
            --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 \
            --INTERLEAVE true -NON_PF true | bwa mem -K 100000000 \
            -v 3 {params.bwa_args} -t {threads} -p  -Y {params.index} /dev/stdin | gatk \
            --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            MergeBamAlignment \
            --VALIDATION_STRINGENCY SILENT \
            --EXPECTED_ORIENTATIONS FR \
            --ATTRIBUTES_TO_RETAIN X0 \
            --ALIGNED_BAM /dev/stdin \
            --UNMAPPED_BAM {input.bam} \
            --OUTPUT {output.bam} \
            -R {params.fa} \
            --PAIRED_RUN true \
            --SORT_ORDER "unsorted" \
            --IS_BISULFITE_SEQUENCE false \
            --ALIGNED_READS_ONLY false \
            --CLIP_ADAPTERS false \
            --MAX_RECORDS_IN_RAM 2000000 \
            --ADD_MATE_CIGAR true \
            --MAX_INSERTIONS_OR_DELETIONS -1 \
            --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
            --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
            --ALIGNER_PROPER_PAIR_FLAGS true \
            --UNMAP_CONTAMINANT_READS true
        '''

# rule dragon_align:
#     input:
#         bam="0_samples/{sample}/{sample}-{unit}.adab.ubam"

#     output:
#         bam="02_alignment-Dragon/{sample}/{sample}-{unit}_mergedUnmapped.bam"

#     conda: "../env/wes_gatk.yml"
#     threads: 4
#     params:
#         fa = ref_fasta,
#         index = ref_bwa,
#         bwa_args = config["aligner_extra_args"]

#     benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
#     resources:
#         mem_mb = 32* 1024,
#         # cores=config["align_threads"],
#         mem_gb = 32,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     shell:
#         '''
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             SamToFastq \
#             -I {input.bam} \
#             --FASTQ /dev/stdout \
#             --CLIPPING_ATTRIBUTE XT --CLIPPING_ACTION 2 \
#             --INTERLEAVE true -NON_PF true | bwa mem -K 100000000 \
#             -v 3 {params.bwa_args} -t {threads} -p  -Y {params.index} /dev/stdin | gatk \
#             --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             MergeBamAlignment \
#             --VALIDATION_STRINGENCY SILENT \
#             --EXPECTED_ORIENTATIONS FR \
#             --ATTRIBUTES_TO_RETAIN X0 \
#             --ALIGNED_BAM /dev/stdin \
#             --UNMAPPED_BAM {input.bam} \
#             --OUTPUT {output.bam} \
#             -R {params.fa} \
#             --PAIRED_RUN true \
#             --SORT_ORDER "unsorted" \
#             --IS_BISULFITE_SEQUENCE false \
#             --ALIGNED_READS_ONLY false \
#             --CLIP_ADAPTERS false \
#             --MAX_RECORDS_IN_RAM 2000000 \
#             --ADD_MATE_CIGAR true \
#             --MAX_INSERTIONS_OR_DELETIONS -1 \
#             --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
#             --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
#             --ALIGNER_PROPER_PAIR_FLAGS true \
#             --UNMAP_CONTAMINANT_READS true
#         '''

rule QC_alignment:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam"

    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/QC_alignment/{sample}.txt"
    output:
        cov = "03_bamPrep/QC/{sample}.cov",
        stats = "03_bamPrep/QC/{sample}.stats"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > {output.cov}
        samtools flagstat {input} > {output.stats}
        """

rule qualimap:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("03_bamPrep/QC/{sample}_Qualimap")
    benchmark: "benchamrks/Qualimap/{sample}.txt"

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (16 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 16  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt  

    shell:
        """
        qualimap \
            bamqc \
            -bam {input} \
            --java-mem-size=15G \
            --paint-chromosome-limits \
            --genome-gc-distr HUMAN \
            -nt {threads} \
            -skip-duplicated \
            --skip-dup-mode 0 \
            -outdir {output} \
            -outformat HTML
        """


## TODO: implement automatic indexing 
# rule index_ref:
#     input: f"{ref_fasta}"
    
#     conda: "../env/wes_gatk.yml"

#     output: f"{ref_fasta}.fai"

#     resources:
#         mem_mb=2048,
#         cores=1,
#         mem_gb=2,
        # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     shell: "samtools faidx {input}"


# rule refrence_dict:
#     input: f"{ref_fasta}"
    
#     conda: "../env/wes_gatk.yml"

#     output: f"{ref_fasta}".replace(".fa", ".dict")
#     resources:
#         mem_mb=2048,
#         cores=1,
#         mem_gb=2,
        # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     shell: "picard CreateSequenceDictionary -R {input}"

# rule bwa_index:
#     input: f"{ref_fasta}"
    
#     conda: "../env/wes_gatk.yml"

#     output: directory(f"{ref_bwa_path}/{ref_prefix}")
#     resources:
#         mem_mb=8192,
#         cores=4,
#         mem_gb=8,
        # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     shell: 
#         "bwa index -p {output} {input}"

## TODO: create input function for align to use either workflows 
# rule bwa_align:
#     input:
#         R1 = "00_trimmomatic/{sample}/{sample}-{unit}1.trimmed.fastq",
#         R2 = "00_trimmomatic/{sample}/{sample}-{unit}2.trimmed.fastq",

    
#     conda: "../env/wes_gatk.yml"

#     output:
#         temp("02_alignment/{sample}/{sample}-{unit}.bam")

#     threads: 4
#     params:
#         fa = ref_fasta,
#         index = ref_bwa 

#         ## TODO: fix this
#         # index = lambda wildcards: ref_bwa if not config["index_fasta"] else rules.bwa_index.output

#         # library_index = lambda wildcards: units.loc[:, 'library_index'][units['unit'] == f"{wildcards.unit}"].tolist()[0],
#         # lane = lambda wildcards: units.loc[:, 'lane'][units['unit'] == f"{wildcards.unit}"].tolist()[0]

#     log: 
#         bwa = "logs/bwa/{sample}/{sample}-{unit}_bwa.log",

#     benchmark: "benchamrks/{sample}/{sample}-{unit}_bwa.txt"
#     resources:
#         mem_mb=int(config["align_mem"])* 1024,
#         cores=config["align_threads"],
#         mem_gb=int(config["align_mem"]),
        # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     ## TODO: why does wildcards have errors in parsing 
#     ##       fix it and use wildcards to generate LB, PL, SM, etc.

#     shell:
#         """
#         R1={input.R1}
#         SM={wildcards.sample}
#         PL="Illumina"
#         LB="{wildcards.sample}_{wildcards.unit}"
#         name=$(basename $R1 | cut -d'_' -f1)
#         RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
#         PU=$RGID.$LB 
#         bwa mem -t {threads} -M \
#             -R "@RG\\tID:$RGID\\tSM:$SM\\tPL:$PL\\tLB:$LB\\tPU:$PU" {params.index} {input.R1} {input.R2} 2> {log.bwa} | \
#             samtools view -1 - >  {output} 
#         """

# rule sort_and_convert_sam:
#     input:
#         "02_alignment/{sample}/{sample}-{unit}.bam"
    
#     conda: "../env/wes_gatk.yml"
#     threads: config["general_low_threads"]
#     output:
#         "02_alignment/{sample}/{sample}-{unit}.sorted.bam"
    
#     resources:
#         mem_mb=int(config["general_low_mem"])* 1024,
#         cores=config["general_low_threads"],
#         mem_gb=int(config["general_low_mem"]),
#         nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         samtools sort {input} -o {output}
#         samtools index {output}
#         """
