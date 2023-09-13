import os
import pandas as pd


#################################################
#                01 define env vars             #
#################################################

# get options
HTDBfile = "/home/marc/Desktop/data/refs/GATK-Resources/Homo_sapiens_assembly38.haplotype_database.txt"
# config["HtDBfile"]

PATH = config["path"]
EXTT = config["ext"]
RS = config["R"]
if config["decompress"]:
    EXT = EXTT.replace(".gz","")
else:
    EXT = EXTT

working_dir = config["working_dir"]
sample_table_file=config.get('sampletable','samples.tsv')
SampleTable = pd.read_table(sample_table_file)

files_R1s = set(SampleTable.iloc[:, 1-1])
files_R2s = set(SampleTable.iloc[:, 8-1])
# sample full names
samples = set(SampleTable.iloc[:, 3-1]) 
units = set(SampleTable.iloc[:, 2-1])
samples_IDs = set(SampleTable.iloc[:, 4-1])
units = (
    pd.read_csv('samples.tsv', sep="\t", dtype={"sample_id": str, "unit": str})
    .set_index(["sample_id", "unit"], drop=False)
    .sort_index()
)

ALL_THREADS = config["threads"]
MEM = config["total_mem"]
GUAP_FOLDER = config["GUAP_DIR"]
R = [1, 2]
source = PATH

working_dir = config["working_dir"]
source_dir = config["GUAP_DIR"]
common_rules = config["common_rules"]

samples_dir = config["input"]

out_dir = config["output"]

ref_bwa = config["reference_index"]

ref_bwa_path = config["reference_output_path"]
ref_prefix = config["reference_output_prefix"]

ref_fasta = config["reference_fasta"]
ref_fasta_path = os.path.dirname(ref_fasta)

bed_file = config["bed_file"]

unique_values = set(line.split('\t')[0] for line in open(bed_file))

known_variants_snps = config["known_variants_snps"]
known_variants_indels = config["known_variants_indels"]
known_variants_indels2 = config["known_variants_indels2"]

nirvana_path = config["nirvana_path"]
annovar_dir = config["annovar_path"]

Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll"
Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/"
Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat"
Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both"

svd_file = config["svd_prefix"]


gff = config["gff_file"]

call_cnv = config["call_CNV"]


#################################################
#            02 define main functions           #
#################################################

def get_raw_gvcf(wildcards):
    return [f"04_calling/RAW/{sample}.gvcf.gz" \
        for sample in set(samples_IDs)]

def get_reblocked_gvcf(wildcards):
    return [f"04-1_gvcf-processing/reblocked/{sample}.gvcf.gz" \
        for sample in set(samples_IDs)]

def get_final_output(wildcards):
    final_output = []

    final_output.extend(expand(
            "04-5_VariantQC/{type}_variants_genotyped_filttered.html",
            type = ["snvs", "indels"]
    ))

    final_output.extend(expand(
            "04_calling/QC/{sample}.eval.grp",
            sample = samples_IDs
    ))

    final_output.extend(expand(
            "04-4_gvcf-QC/{sample}.variant_calling_summary_metrics",
            sample = samples_IDs
    ))
    
    final_output.extend(expand(
            "03-3_bamQC/ValidateSamFile/{sample}.txt",
            sample = samples_IDs
    ))

    final_output.extend(expand(
            "03-3_bamQC/SequencingArtifactToOxoG/{sample}.oxog_metrics",
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
            "05_Annotation/snpEff/{type}_annotation.csv",
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
        final_output.extend(expand(
            "08_cnv_postprocessing/{sample}.intervals_cohort.vcf.gz",
            sample = samples_IDs
    ))
    return final_output

def get_merge_input(wildcards):
    return expand(
            "02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam",
            unit=units.loc[wildcards.sample, "unit"].tolist(),
            sample=wildcards.sample
        )

# for rule all!
rule multiqc:
    input:
        get_final_output
    
    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/Multiqc/report.txt"

    output:
        "multiqc/multiqc_report.html"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        "multiqc . -o multiqc/"

#################################################
#            03 sample pre-processing           #
#################################################

rule FastqToSam:
    input:
        R1 = f"{samples_dir}/{{sample}}-{{unit}}1.{EXT}",
        R2 = f"{samples_dir}/{{sample}}-{{unit}}2.{EXT}"
    
    conda: "../env/wes_gatk.yml"

    output:
        ubam = "0_samples/{sample}/{sample}-{unit}.ubam"
    threads: 4
    params:
        ext = EXT
    benchmark: "benchamrks/FastqToSam/{sample}/{sample}-{unit}.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt

    shell:
        """
        ext={params.ext}
        R1={input.R1}
        SM={wildcards.sample}
        PL="Illumina"
        LB=$SM
        if [[ "$ext" == *".gz" ]]; then
            RGID=$(head -n1 <(zcat {input.R1}) | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
            #RGID=$(zcat {input.R1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        else
            RGID=$(head {input.R1} -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        fi
        ## TODO: confirm to use this
        PU=$RGID.$LB 
        
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            FastqToSam \
            -F1 {input.R1} \
            -F2 {input.R2} \
            -O {output} \
            -SM $SM \
            -LB $LB \
            -PL $PL \
            -RG $RGID \
            -PU $PU \
        """

rule MarkIlluminaAdapters:
    input:
        "0_samples/{sample}/{sample}-{unit}.ubam"
    output:
        bam="0_samples/{sample}/{sample}-{unit}.adab.ubam",
        metrics="0_samples/{sample}/{sample}-{unit}.adap_metrics.txt"
    threads:1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt

    benchmark: "benchamrks/MarkIlluminaAdapters/{sample}/{sample}-{unit}.txt"

    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            MarkIlluminaAdapters \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            # > {log} 2>&1
        '''

#################################################
#            04 mapping to reference            #
#################################################
"""
--ATTRIBUTES_TO_REMOVE NM:
The --ATTRIBUTES_TO_REMOVE parameter specifies a list of attributes to remove from the read records in the BAM file. In this case, NM refers to the "Edit distance to the reference" attribute, which indicates the number of mismatches and gaps in the alignment between a read and the reference genome. Removing this attribute could impact downstream analyses that rely on this information, such as variant calling and quality assessment.

Implication: Removing the NM attribute might affect the quality of variant calling and other analyses that depend on alignment information. It's generally advisable to retain this attribute unless you have a specific reason to remove it.

--ATTRIBUTES_TO_REMOVE MD:
The MD attribute stands for "Mismatching positions," which indicates the positions in the read that differ from the reference sequence. Similar to the NM attribute, removing the MD attribute could impact variant calling and alignment-related analyses.

Implication: Removing the MD attribute may affect alignment-based analyses. Keeping this attribute is usually recommended.

--ADD_PG_TAG_TO_READS false:
The --ADD_PG_TAG_TO_READS parameter controls whether Program Group (PG) tags should be added to the read groups of the output BAM file. PG tags provide information about the program and version used to generate the data, which can be useful for tracking the processing history of the data.

Implication: If you set this parameter to false, PG tags will not be added to the read groups in the output BAM file. This might not have a significant impact on the variant calling process itself, but it could affect downstream data analysis workflows that rely on accurate metadata tracking.

MergeBamAlignment step ensures that the attributes from both sets of reads are properly integrated. The purpose of removing these attributes could be related to consistency and avoiding redundancy, as merging may introduce duplicates or conflicts in attribute values. However, removing them could affect downstream analyses that rely on accurate alignment and mismatch information.

Regarding whether you need to rerun the analysis if you didn't use these parameters: It depends on the specific goals of your analysis and the downstream analyses you plan to perform. If you have already performed variant calling without these parameters and are satisfied with the results, you may not necessarily need to rerun the analysis. However, it's a good practice to carefully consider the implications of removing important attributes and metadata from the BAM file, as it could affect the reliability and interpretability of your results.

When using the Haplotype Caller (part of GATK), the impact of these parameters could be similar, as the principles of alignment and variant calling still apply. However, the specific implications might vary depending on the algorithm's behavior and the downstream analyses you plan to perform on the haplotype-called variants. It's always a good idea to refer to the GATK documentation and best practices for guidance on parameter choices and their effects on different analysis scenarios.
"""

rule ubam_align:
    input:
        bam="0_samples/{sample}/{sample}-{unit}.adab.ubam"

    output:
        bam="02_alignment/{sample}/{sample}-{unit}_mergedUnmapped.bam"

    conda: "../env/wes_gatk.yml"
    threads:8
    params:
        fa = ref_fasta,
        index = ref_bwa,
        bwa_args = config["aligner_extra_args"],
        warp_workflow_params = "--ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD --ADD_PG_TAG_TO_READS false"

    benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
    resources:
        mem_mb = 24* 1024,
        mem_gb = 24,
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
            --UNMAP_CONTAMINANT_READS true {params.warp_workflow_params}
        '''

#################################################
#                05 mapping QC                  #
#################################################


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
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 2  * attempt,
        runtime = lambda wildcards, attempt: 60 * 1 * attempt
    shell:
        """
        samtools depth {input} | awk '{{sum+=$3}} END {{print "Average = ",sum/NR, "No of covered Nuc = ", NR}}' > {output.cov}
        samtools flagstat {input} > {output.stats}
        """

rule qualimap:
    input:
        bam = "03_bamPrep/{sample}.dedub.sorted.bam",
        bed = "resources/padded.bed"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("03_bamPrep/QC/{sample}_Qualimap")
    benchmark: "benchamrks/Qualimap/{sample}.txt"

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (10 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 10  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt  

    shell:
        """
        qualimap \
            bamqc \
            -bam {input.bam} \
            --java-mem-size=15G \
            --paint-chromosome-limits \
            --genome-gc-distr HUMAN \
            -nt {threads} \
            -skip-duplicated \
            --skip-dup-mode 0 \
            -outdir {output} \
            --feature-file {input.bed} \
            -outformat HTML
        """

#################################################
#              06 bam processing                #
#################################################

rule MergeSamFiles:
    input:
        get_merge_input
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Mergesamfiles/{sample}.txt"

    output:
        "03_bamPrep/{sample}.bam"
    params:
        bams = lambda wildcards: [f" -I 02_alignment/{wildcards.sample}/{wildcards.sample}-{b}_mergedUnmapped.bam" for b in units.loc[wildcards.sample, "unit"].tolist()],
        ref = ref_fasta

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt
    shell:
        """
        picard MergeSamFiles  \
            {params.bams} \
            -OUTPUT {output} \
            --SORT_ORDER "unsorted" 
        """


rule MarkDuplicates:
    input:
        "03_bamPrep/{sample}.bam"
    
    conda: "../env/wes_gatk.yml"

    output:
        bam = "03_bamPrep/{sample}.dedub.bam",
        matrix = "03_bamPrep/{sample}.dedub.matrix"

    threads: 1
    benchmark: "benchamrks/mrkDuplicates/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        picard MarkDuplicates -I {input} \
            -O {output.bam} \
            -M {output.matrix} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true \
            --CLEAR_DT false \
            --ADD_PG_TAG_TO_READS false
        """

rule SortNFix:
    input:
        "03_bamPrep/{sample}.dedub.bam"
    conda: "../env/wes_gatk.yml"
    output:
        "03_bamPrep/{sample}.dedub.sorted.bam"

    threads: 4
    params:
        fa = ref_fasta

    benchmark: "benchamrks/SortNFix/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (10 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 10  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
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

rule BaseRecalibrator:
    input: "03_bamPrep/{sample}.dedub.sorted.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}.report"
    params: 
        known_sites = known_variants_snps,
        known_sites2 = known_variants_indels,
        known_sites3 = known_variants_indels2,
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]

    threads: 1
    benchmark: "benchamrks/BaseRecalibrator/{sample}.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            BaseRecalibrator \
            -R {params.ref} \
            -I {input} \
            --use-original-qualities \
            --known-sites {params.known_sites} \
            --known-sites {params.known_sites2} \
            --known-sites {params.known_sites3} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -O {output} 
        """

rule applyBaseRecalibrator:
    input: 
        bam = "03_bamPrep/{sample}.dedub.sorted.bam",
        report = "03_bamPrep/{sample}.report"
    benchmark: "benchamrks/ApplyBaseRecab/{sample}.txt"
    
    conda: "../env/wes_gatk.yml"

    output: "03_bamPrep/{sample}.pqsr.bam"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    params: 
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ApplyBQSR  -R {params.ref} \
            -I {input.bam} \
            -bqsr {input.report} -O {output} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-md5 \
            --use-original-qualities

        samtools index {output} 
        """

rule bqsr_calibrated_report:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/recaliprated_report/{sample}.txt"

    output: "03_bamPrep/{sample}_pqsr.report"
    params: 
        known_sites = known_variants_snps,
        known_sites2 = known_variants_indels,
        known_sites3 = known_variants_indels2,
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"]
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" BaseRecalibrator -R {params.ref} \
            -I {input} \
            --known-sites {params.known_sites} \
            --known-sites {params.known_sites2} \
            --known-sites {params.known_sites3} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -O {output} 
        """

rule AnalyzeCovariates:
    input: 
        raw = "03_bamPrep/{sample}.report", 
        bqsr = "03_bamPrep/{sample}_pqsr.report"

    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/analyzeCovariates/{sample}.txt"

    output: "03_bamPrep/QC/{sample}.pdf"

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" AnalyzeCovariates -before {input.raw} \
            -after {input.bqsr} -plots {output}
        """

#-----------------------------------------------#
#               get contamination               #
#-----------------------------------------------#

rule create_padded_interval:
    input: 
        bed = bed_file

    conda: "../env/wes_gatk.yml"
    params:
        padding = config["padding"],
        ref = ref_fasta

    output: "resources/padded.bed"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 2  * attempt,
        runtime = lambda wildcards, attempt: 30 * attempt

    shell:
        """
        bedtools slop -i {input.bed} -g {params.ref}.fai -b {params.padding} > {output}
        """

rule GenerateSubsettedContaminationResources:
    input: 
        intervalslist = "resources/padded.bed",
        UD = svd_file + ".UD",
        BED = svd_file + ".bed",
        MU = svd_file + ".mu",

    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/GenSubCont/GenerateSubsettedContaminationResources.txt"

    output: 
        UD = "03-1_SubsettedContamination/EXOME_Contams.UD",
        BED = "03-1_SubsettedContamination/EXOME_Contams.bed",
        MU = "03-1_SubsettedContamination/EXOME_Contams.mu",
        target_overlap_counts = "03-1_SubsettedContamination/target_overlap_counts.txt"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 2  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt

    shell:
        '''
        grep -vE "^@" {input.intervalslist} |
            awk -v OFS='\t' '$2=$2-1' |
            bedtools intersect -c -a {input.BED} -b - |
            cut -f6 > {output.target_overlap_counts}

        function restrict_to_overlaps() {{
            # print lines from whole-genome file from loci with non-zero overlap
            # with target intervals
            WGS_FILE=$1
            EXOME_FILE=$2
            paste {output.target_overlap_counts} $WGS_FILE |
                grep -Ev "^0" |
                cut -f 2- > $EXOME_FILE
            echo "Generated $EXOME_FILE"
        }}

        restrict_to_overlaps {input.UD} {output.UD}
        restrict_to_overlaps {input.BED} {output.BED}
        restrict_to_overlaps {input.MU} {output.MU}

        '''


rule CheckContamination:
    input: 
        UD = "03-1_SubsettedContamination/EXOME_Contams.UD",
        BED = "03-1_SubsettedContamination/EXOME_Contams.bed",
        MU = "03-1_SubsettedContamination/EXOME_Contams.mu",
        target_bam = "03_bamPrep/{sample}.bam",

        ref = ref_fasta,

    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/GenSubCont/CheckContamination.{sample}.txt"

    output: 
        selfSM = "03-2_contamitnation_check/{sample}.selfSM",
        contamination = "03-2_contamitnation_check/{sample}.contamination.txt"  
    params:
        get_con_file = f"{source_dir}/workflow/scripts/get_contamination.py",
        prefix = "03-2_contamitnation_check/{sample}",
        svdprefix = "03-1_SubsettedContamination/EXOME_Contams"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (2 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 2  * attempt,
        runtime = lambda wildcards, attempt: 30 * attempt

    shell:
        '''
        python3 {params.get_con_file} \
            --prefix {params.prefix} \
            --output {output.contamination} \
            --svdprefix {params.svdprefix} \
            --get_con_file {params.get_con_file} \
            --target_bam {input.target_bam} \
            --ref {input.ref} \
            --threads {threads} 
        '''


#################################################
#             07 variant callign                #
#################################################

# TODO:  update once defined this rule!
Dragon_STR_Model_path = ""

rule HaplotypeCaller:
    input: 
        bam = "03_bamPrep/{sample}.pqsr.bam",
        contamination = "03-2_contamitnation_check/{sample}.contamination.txt"
    
    conda: "../env/wes_gatk.yml"

    output: 
        vcf = "04_calling/RAW/{sample}.gvcf.gz",
        bam = "04_calling/BAMS/{sample}_raw.gvcf.bamout.bam",
    
    params: 
        ref = ref_fasta,
        bed = bed_file,
        extra_args = config["caller_extra_args"],
        padding = config["padding"],
        Dragon_mode = lambda wildcards: f"--dragen-mode --disable-spanning-event-genotyping" if config['Dragon_mode'] else "",
        Dragon_STR = lambda wildcards: f"--dragstr-params-path {Dragon_STR_path}" if config['Dragon_mode'] else "",

    benchmark: "benchamrks/HaplotypeCaller/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (12 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 12  * attempt,
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
        reduced = lambda wildcards, attempt: int(attempt * (12 * 0.80))

    threads: 1

    shell:
        """
        cont=$(cat {input.contamination})
        echo $cont
        
        gatk --java-options "-Xmx{resources.reduced}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller -R {params.ref} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation {params.extra_args} \
            -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -I {input.bam} --native-pair-hmm-threads {threads} -ERC GVCF -O {output.vcf} \
            -contamination $cont \
            -bamout {output.bam} \

        """

rule ReblockGVCF:
    input: "04_calling/RAW/{sample}.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04-1_gvcf-processing/reblocked/{sample}.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        padding = config["padding"],
        tree_score_threashold = "--tree-score-threshold-to-no-call N"
        # TODO: update the value and import it in the pipeline shell 


    benchmark: "benchamrks/Reblock/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        

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
            --floor-blocks -OVI \
            -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -do-qual-approx \
            -V {input} -O {output}

        """

#-----------------------------------------------#
#           DRAGON-specific rules               #
#-----------------------------------------------#

# which turned out to be filteration for single 
# sample runs as said here:
# https://www.melbournebioinformatics.org.au/tutorials/tutorials/variant_calling_gatk1/variant_calling_gatk1/

rule HardFilterVcf:
    input: "04-1_gvcf-processing/reblocked/{sample}.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04-1_gvcf-processing/filtterd/{sample}_reblocked.gvcf.gz"

    params: 
        bed = bed_file,
        padding = config["padding"],
        Dragon_filter = "--filter-expression 'QUAL < 10.4139' --filter-name 'DRAGENHardQUAL' "


    benchmark: "benchamrks/Reblock/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            VariantFiltration \
            -V {input} \
            -L {params.bed} \
            --interval-padding {params.padding} \
            --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
            --filter-name "HardFiltered" \
            -O {output}
        """

rule CNNScoreVariants:
    input: 
        vcf = "04-1_gvcf-processing/filtterd/{sample}_reblocked.gvcf.gz",
        bam = "04_calling/BAMS/{sample}_raw.gvcf.bamout.bam",
    
    conda: "../env/wes_gatk.yml"

    output: "04-1_gvcf-processing/CNNScoreVariants/{sample}.gvcf.gz"

    params: 
        bed = bed_file,
        padding = config["padding"],
        ref = ref_fasta,

    benchmark: "benchamrks/Reblock/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CNNScoreVariants \
            -V {input.vcf} \
            -R {params.ref} \
            -O {output} \
            -I {input.bam} \
            -tensor-type read-tensor
        """

rule FilterVariantTranches:
    input: 
        "04-1_gvcf-processing/CNNScoreVariants/{sample}.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: 
        "04-1_gvcf-processing/{sample}_filterred.gvcf.gz"

    params: 
        bed = bed_file,
        padding = config["padding"],
        ref = ref_fasta,
        hapmap_resource_vcf = "hapmap_3.3.hg38.vcf.gz",
        omni_resource_vcf = "1000G_omni2.5.hg38.vcf.gz",
        one_thousand_genomes_resource_vcf = "1000G_phase1.snps.high_confidence.hg38.vcf.gz",
        dbsnp_resource_vcf = "Homo_sapiens_assembly38.dbsnp138.vcf",
        # info_key = "",
        # snp_tranches = {sep=" " prefix("--snp-tranche ", snp_tranches)} 
        # indel_tranches = {sep=" " prefix("--indel-tranche ", indel_tranches)} 
        
    benchmark: "benchamrks/Reblock/{sample}.txt"

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            FilterVariantTranches \
            -V ~{input_vcf} \
            -O ~{vcf_basename}.filtered.vcf.gz \
            --resource {params.hapmap_resource_vcf} \
            --resource {params.omni_resource_vcf} \
            --resource {params.one_thousand_genomes_resource_vcf} \
            --resource {params.dbsnp_resource_vcf} \
            --info-key CNN_1D \
            --snp-tranche 99.95 \
            --indel-tranche 99.4 \
            --create-output-variant-index true
        """

#################################################
#                08 Genotyping                  #
#################################################

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
            echo -e  "${{i}}\t04-1_gvcf-processing/reblocked/${{i}}.gvcf.gz" | cat >> {output}
        done
        """

rule create_new_interval:
    input: 
        bedfile = bed_file
    output:
        newbed = "resources/Processed_intervals.bed"

    threads: 1

    conda: "../env/wes_gatk.yml"

    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
        

    shell:
        """
        sort -k1,1 -k2,2n {input.bedfile} | awk '
        {{
            # If it is a different chromosome or if it is a non-overlapping interval
            if ($1 != chrom) {{
                # If it is not the first interval
                if (chrom) {{
                    print chrom, start, end;
                }}
                chrom = $1;
                start = $2;
                end = $3;
            }} else if ($3 > end) {{
                # If it is an overlapping interval
                end = $3;
            }}
        }}
        END {{
            # Print the last interval
            print chrom, start, end;
        }}
        ' > {output.newbed}
        """

rule scatter_intervals:
    input:
        ref = ref_fasta,
        intervals = "resources/Processed_intervals.bed"
    
    conda: "../env/wes_gatk.yml"

    output:
        dir = directory("resources/Scatterred_intervals"),
        intervals = expand("resources/Scatterred_intervals/00{interval}-scattered.interval_list", interval = [f"{i:02d}" for i in range(0, int(config['scatter_count']))])

    params:
        scatter_count = int(config['scatter_count']),

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
        
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            SplitIntervals \
            -R {input.ref} \
            -L {input.intervals} \
            --scatter-count {params.scatter_count} \
            --interval-merging-rule ALL \
            -mode INTERVAL_SUBDIVISION \
            -O {output.dir} 
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
        runtime = lambda wildcards, attempt: 60 * 4 * attempt,
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
            --genomicsdb-shared-posixfs-optimizations true \

        """

rule genotype_gvcfs:
    input:
        infile = "04-2_GenomicsDB/Interval_{interval}",
        intervals = "resources/Scatterred_intervals/00{interval}-scattered.interval_list"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

    output:
        "04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz"
    params:
        ref = ref_fasta,
        Additional_annotation = "",

    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            GenotypeGVCFs \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            -R {params.ref} -V gendb://{input.infile} -O {output} \
            {params.Additional_annotation} --merge-input-intervals \
            -L {input.intervals} \
            --only-output-calls-starting-in-intervals \

        """

rule combine_gvcf:
    input:
        gvcfs= expand("04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz", interval = [f"{i:02d}" for i in range(0, int(config['scatter_count']))])
    
    conda: "../env/wes_gatk.yml"

    output:
        "04_calling/variants_genotyped.gvcf.gz"

    params:
        gvcfs = lambda wildcards, input: [f" --input {v}" for v in input["gvcfs"]],
        ref = ref_fasta

    benchmark: "benchamrks/combine_gvcf/combinedvcf.txt"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
            GatherVcfsCloud \
            --ignore-safety-checks \
            --gather-type BLOCK \
            {params.gvcfs} \
            -O {output}

        tabix {output}
        """

#################################################
#            09 Variant processing              #
#################################################

rule VariantEval:
    input: "04_calling/RAW/{sample}.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/QC/{sample}.eval.grp"
    threads: 1
    params:
        ref = ref_fasta
    benchmark: "benchamrks/VariantEval/{sample}/rawvcf.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        touch {output}
        gatk VariantEval \
            -R {params.ref} \
            --eval:{wildcards.sample} {input} \
            -O {output}
        """

rule Split_variants_idnel:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Split_variants_idnel/genotypedvcf.txt"

    output: "04_calling/indels/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include INDEL \
            -O {output}
        """

rule variant_filteration_indels:
    input: 
        "04_calling/indels/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/variant_filteration_indels/indels.txt"

    output: 
        "04_calling/indels/variants_genotyped.filttered.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk VariantFiltration \
            --variant {input} \
            --filter-expression "QD < 2.0"                  --filter-name "QD2" \
            --filter-expression "SOR > 10.0"                  --filter-name "SOR10" \
            --filter-expression "FS > 200.0"                --filter-name "FS200" \
            --filter-expression "ReadPosRankSum < -20.0"    --filter-name "ReadPosRankSum-20" \
            --filter-expression "InbreedingCoeff < -0.8"    --filter-name "InbreedingCoeff-0.8" \
            --create-output-variant-index true \
            --output {output}
        """


rule Split_variants_snp:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Split_variants_snp/genotypedvcf.txt"

    output: "04_calling/snvs/variants_genotyped.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}"  SelectVariants \
            -R {params.ref} \
            -V {input} \
            --select-type-to-include SNP \
            -O {output}
        """

rule variant_filteration_snps:
    input: "04_calling/snvs/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/variant_filteration_snps/snps.txt"

    output: "04_calling/snvs/variants_genotyped.filttered.gvcf.gz"
    params:
        ref = ref_fasta
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 2
    shell:
        """
        gatk VariantFiltration \
            --variant {input} \
            --filter-expression "QD < 2.0"              --filter-name "QD2" \
            --filter-expression "SOR > 3.0"             --filter-name "SOR3" \
            --filter-expression "FS > 60.0"             --filter-name "FS60" \
            --filter-expression "MQ < 40.0"             --filter-name "MQ40" \
            --filter-expression "MQRankSum < -12.5"     --filter-name "MQRankSum-12.5" \
            --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
            --create-output-variant-index true \
            --output {output}
        """

#################################################
#            10 Variant annotation              #
#################################################

rule consequence:
    input: "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/consequence/{type}/bcftools_csq.txt"

    output: "04_calling/QC/{type}/bcftools_csq.vcf"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads:1
    params:
        ref = ref_fasta,
        gff = gff

    shell:
        """
        bcftools csq -f {params.ref} -g {params.gff} {input} -Ov -o {output}
        """

rule bcftools_stats:
    input: "04_calling/QC/{type}/bcftools_csq.vcf"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/bcftools_stats/{type}/bcftools_csq.txt"

    output: "04_calling/QC/{type}/bcftools.stats"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
            bcftools stats {input} > {output}
        """

rule plot_bcftools_stats:
    input: "04_calling/QC/{type}/bcftools.stats"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/plot_bcftools_stats/{type}/bcftools_plots.txt"
    output: directory("04_calling/QC/{type}/bcftools_plots")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        plot-vcfstats -p {output} {input}
        """

rule Nirvana:
    input:
        "04_calling/{type}/variants_genotyped.filttered.gvcf.gz"
    
    benchmark: "benchamrks/Nirvana/{type}/variants_genotyped_annotation.txt"

    output:
        "05_Annotation/Nirvana/{type}/Annotation.json.gz"
    threads: 1

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    params:
        Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/",
        Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat",
        Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both",
        Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll",
        file_name = lambda wildcards: f"05_Annotation/Nirvana/{wildcards.type}/Annotation",
        extra = config["nirvana_extra_args"]

    shell:
        """
        eval "$(conda shell.bash hook)"
        set +u; conda activate dotnet; set -u
        dotnet {params.Nirvana_cmd} \
            -i {input} \
            -o {params.file_name} \
            -r {params.Nirvana_ref} \
            --sd {params.Nirvana_supplementray} \
            -c {params.Nirvana_cache} {params.extra}
        """

rule snpEFF:
    input:
        "04_calling/{type}/variants_genotyped.filttered.gvcf.gz"

    conda: "../env/wes_gatk.yml"   
    benchmark: "benchamrks/snpeff/{type}/variants_genotyped_annotation.txt"

    output:
        stats = "05_Annotation/snpEff/{type}_annotation.csv",
        html = "05_Annotation/snpEff/{type}_annotation.html",
        vcf = "05_Annotation/snpEff/{type}_annotation.vcf"

    threads: 1

    resources:
        mem_mb=lambda wildcards, attempt: (32 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 32  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    params:
        genome =  config["snpeff_genome"]

    shell:
        """
        snpEff -Xmx{resources.mem_gb}G \
            -v -csvStats {output.stats} \
            -stats {output.html} \
            {params.genome} {input} > {output.vcf}
        """

rule Annovar:
    input:
        "04_calling/{type}/variants_genotyped.filttered.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Annovar/{type}/variants_genotyped_annotation.txt"
    output:
        directory("05_Annotation/ANNOVAR/{type}")
    threads: 32
    resources:
        mem_mb=32 * 1024,
        mem_gb=32,
        runtime = lambda wildcards, attempt: 60 * 24 * attempt
    params:
        annovar_dir = annovar_dir,
        protocol = config["annovar_protocol"],
        operation = config["annovar_operation"],
        output = lambda wildcards: f"05_Annotation/ANNOVAR/{wildcards.type}/annotations",
        extra = config["annovar_extra_args"]

    shell:
        """
        mkdir -p {output}
        perl {params.annovar_dir}/table_annovar.pl {input} {params.annovar_dir}/humandb/ \
            -buildver hg38 \
            -out {params.output} -remove \
            -protocol {params.protocol} \
            -operation {params.operation} \
            -nastring . \
            -vcfinput \
            --thread {threads} {params.extra}
        """

#-----------------------------------------------#
#                    BAM QC                     #
#-----------------------------------------------#

# Input is the final processed bam file

# Collect sequencing yield quality metrics
rule CollecQualityYieldMetrics:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam",
    output:
        "03-3_bamQC/{sample}.metrics.txt",

    conda: "../env/wes_gatk.yml"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" 
            CollectQualityYieldMetrics \
            --INPUT {input} \
            --OQ true \
            --OUTPUT {output}
        """

# Collect alignment summary and GC bias quality metrics
rule CollectBamQualityMetrics:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam",
    output:
        pre_adapter_detail_metrics = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix.pre_adapter_detail_metrics",
        pre_adapter_summary_metrics = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix.pre_adapter_summary_metrics",

    conda: "../env/wes_gatk.yml"
    params: 
        prefix = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix",
        ref = ref_fasta,

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CollectMultipleMetrics \
            --INPUT {input} \
            --REFERENCE_SEQUENCE {params.ref} \
            --OUTPUT {params.prefix} \
            --ASSUME_SORTED true \
            --PROGRAM null \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectSequencingArtifactMetrics \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM CollectGcBiasMetrics \
            --METRIC_ACCUMULATION_LEVEL null \
            --METRIC_ACCUMULATION_LEVEL SAMPLE \
            --METRIC_ACCUMULATION_LEVEL LIBRARY \
            --METRIC_ACCUMULATION_LEVEL ALL_READS 

        """



# Collect ConvertSequencingArtifactToOxoG
rule ConvertSequencingArtifactToOxoG:
    input:
        bam = "03_bamPrep/{sample}.dedub.sorted.bam",
        pre_adapter_detail_metrics = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix.pre_adapter_detail_metrics",
    output:
        "03-3_bamQC/SequencingArtifactToOxoG/{sample}.oxog_metrics"

    conda: "../env/wes_gatk.yml"
    params: 
        prefix = "03-3_bamQC/SequencingArtifactToOxoG/{sample}",
        inPasename = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix",
        ref = ref_fasta,
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ConvertSequencingArtifactToOxoG \
            --INPUT_BASE {params.inPasename} \
            --OUTPUT_BASE {params.prefix} \
            --REFERENCE_SEQUENCE {params.ref} \
        """


# Check that the fingerprints of separate readgroups all match
rule CrossCheckFingerprints:
    input:
        bams = expand("03_bamPrep/{sample}.dedub.sorted.bam", sample = samples_IDs),

    output:
        "03-3_bamQC/CrossCheckFingerprints/metrics.file"

    conda: "../env/wes_gatk.yml"
    params: 
        files = lambda wildcards, input: [f" --INPUT {v}" for v in input["bams"]],
        ref = ref_fasta,
        htdbf = HTDBfile
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CrosscheckFingerprints \
            --OUTPUT {output} \
            --HAPLOTYPE_MAP {params.htdbf} \
            --EXPECT_ALL_GROUPS_TO_MATCH true \
            {params.files} \
            --LOD_THRESHOLD -10.0 \
        """


rule CollectQualityYieldMetrics:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam",
    output:
        directory("03-3_bamQC/CollectQualityYieldMetrics/{sample}/"),

    conda: "../env/wes_gatk.yml"
    params: 
        prefix = "03-3_bamQC/CollectQualityYieldMetrics/{sample}/BamMatrix",
        ref = ref_fasta,

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CollectMultipleMetrics \
            --INPUT {input} \
            --REFERENCE_SEQUENCE {params.ref} \
            --OUTPUT {params.prefix} \
            --ASSUME_SORTED true \
            --PROGRAM null \
            --PROGRAM CollectBaseDistributionByCycle \
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM MeanQualityByCycle \
            --PROGRAM QualityScoreDistribution \
            --METRIC_ACCUMULATION_LEVEL null \
            --METRIC_ACCUMULATION_LEVEL ALL_READS

        """

rule CollectReadgroupBamQualityMetrics:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam",
    output:
        directory("03-3_bamQC/CollectReadgroupBamQualityMetrics/{sample}/"),

    conda: "../env/wes_gatk.yml"
    params: 
        prefix = "03-3_bamQC/CollectReadgroupBamQualityMetrics/{sample}/RGBamMatrix",
        ref = ref_fasta,

    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CollectMultipleMetrics \
            --INPUT {input} \
            --REFERENCE_SEQUENCE {params.ref} \
            --OUTPUT {params.prefix} \
            --ASSUME_SORTED true \
            --PROGRAM null \
            --PROGRAM CollectAlignmentSummaryMetrics \
            --PROGRAM CollectGcBiasMetrics \
            --METRIC_ACCUMULATION_LEVEL null \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP
        """

rule ValidateSamFile:
    input:
        "03_bamPrep/{sample}.dedub.sorted.bam",
    output:
        "03-3_bamQC/ValidateSamFile/{sample}.txt",

    conda: "../env/wes_gatk.yml"
    params: 
        ref = ref_fasta

    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ValidateSamFile \
            --INPUT {input} \
            --OUTPUT {output} \
            --REFERENCE_SEQUENCE {params.ref} \
            --MODE VERBOSE \
            --IS_BISULFITE_SEQUENCED false
        """

#-----------------------------------------------#
#                    VCF QC                     #
#-----------------------------------------------#

rule ValidateVariants:
    input:
        "04-1_gvcf-processing/reblocked/{sample}.gvcf.gz",
    conda: "../env/wes_gatk.yml"
    params: 
        ref = ref_fasta

    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            ValidateSamFile \
            -V {input} \
            -R {params.ref} \
            --validation-type-to-exclude ALLELES
        """ 

rule BedToIntervalList:
    input:
        bed = "resources/padded.bed",
        seq_dict = str(config["reference_fasta"]).replace(".fasta", ".dict")
    output:
        "resources/padded.intervals"
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    conda: "../env/wes_gatk.yml"

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            BedToIntervalList \
            -I {input.bed} \
            -O {output} \
            -SD {input.seq_dict}
        """

rule CollectVariantCallingMetrics:
    input:
        vcf = "04-1_gvcf-processing/reblocked/{sample}.gvcf.gz",
        intervals = "resources/padded.intervals"
    output:
        summary = "04-4_gvcf-QC/{sample}.variant_calling_summary_metrics",
        detail_metrics = "04-4_gvcf-QC/{sample}.variant_calling_detail_metrics",
    params: 
        prefix = "04-4_gvcf-QC/{sample}",
        refdict = ref_fasta.replace(".fasta", ".dict"),
        dbsnp = known_variants_snps,
        
    
    resources:
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
    threads: 1
    conda: "../env/wes_gatk.yml"

    shell:
        """
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            CollectVariantCallingMetrics \
            -I {input.vcf} \
            --OUTPUT {params.prefix} \
            --DBSNP {params.dbsnp} \
            --SEQUENCE_DICTIONARY {params.refdict} \
            --TARGET_INTERVALS {input.intervals} \
            --GVCF_INPUT true
        """ 

rule VariantQC:
    input:
        "04_calling/{type}/variants_genotyped.filttered.gvcf.gz",
    output:
        "04-5_VariantQC/{type}_variants_genotyped_filttered.html",
    params: 
#        path_to_tool = "/mnt/home/mansourt/annDB/DISCVRSeq",
        path_to_tool = "/mnt/home/mansourt/annDB",
        ref = ref_fasta
        
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 24 * attempt
    threads: 1
    conda: "../env/wes_gatk.yml"

    shell:
        """
        java \
            -jar {params.path_to_tool}/DISCVRSeq-1.3.49.jar VariantQC \
            -R {params.ref} \
            -V {input} \
            -O {output}
        """

#--------------------------------------------------------------------------------------------------------------------------#
#                                       E N D   O F   M A I N   P i p e L i n e                                            #
#--------------------------------------------------------------------------------------------------------------------------#



# #-----------------------------------------------#
# #               Dragon tasks                    #
# #-----------------------------------------------#
# rule dragon_ref_prep:
#     input:
#         ref = ref_fasta,
#     output:
#         directory("resources/dragon_ref/hg38_no_alt")

#     conda: "../env/wes_gatk.yml"
#     threads: 4

#     benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
#     resources:
#         mem_mb = 32* 1024,
#         mem_gb = 32,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         '''
#         outdir=resources/dragon_ref/hg38_no_alt
#         mkdir -p $outdir
#         dragen-os --build-hash-table true --ht-reference {input}  --output-directory $outdir --num-threads {threads}
#         '''

# rule dragon_align:
#     input:
#         bam="0_samples/{sample}/{sample}-{unit}.adab.ubam",
#         ref = ref_fasta,
#         dragon_ref = "resources/dragon_ref/hg38_no_alt"

#     output:
#         bam="02_alignment/{sample}/{sample}-{unit}_mergedDragonUnmapped.bam"

#     conda: "../env/wes_gatk.yml"
#     threads: 4
#     params:
#         fa = ref_fasta,
#         index = ref_bwa,
#         bwa_args = config["aligner_extra_args"],
#         daragon_workflow_params = "--ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD --ADD_PG_TAG_TO_READS false"

#     benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
#     resources:
#         mem_mb = 32* 1024,
#         mem_gb = 32,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     shell:
#         '''
#         dragen-os \
#             -b {input.bam} \
#             -r {input.dragon_ref} --num-threads {threads} \
#             --interleaved=1 | samtools view \
#             -h -O BAM - > aligned.bam
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             MergeBamAlignment \
#             --VALIDATION_STRINGENCY SILENT \
#             --EXPECTED_ORIENTATIONS FR \
#             --ATTRIBUTES_TO_RETAIN X0 \
#             --ATTRIBUTES_TO_REMOVE RG \
#             --ATTRIBUTES_TO_REMOVE NM \
#             --ATTRIBUTES_TO_REMOVE MD \
#             --ALIGNED_BAM aligned.bam \
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
#             --UNMAP_CONTAMINANT_READS true {params.daragon_workflow_params}
#             --ADD_PG_TAG_TO_READS false
#         rm aligned.bam

#         '''


#-----------------------------------------------#
#                   CNVs                        #
#-----------------------------------------------#


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
