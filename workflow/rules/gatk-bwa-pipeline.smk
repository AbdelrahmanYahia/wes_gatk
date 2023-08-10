import os
import pandas as pd


#################################################
#                01 define env vars             #
#################################################

# get options
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
            "04_calling/QC/{sample}.eval.grp",
            sample = samples_IDs
    ))
    final_output.extend(expand(
            "03-2_contamitnation_check/{sample}.contamination.txt",
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
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

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
    threads:4
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

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
        bwa_args = config["aligner_extra_args"],
        daragon_workflow_params = "--ATTRIBUTES_TO_REMOVE NM --ATTRIBUTES_TO_REMOVE MD --ADD_PG_TAG_TO_READS false"

    benchmark: "benchamrks/ubam_align/{sample}/{sample}-{unit}.txt"
    resources:
        mem_mb = 32* 1024,
        # cores=config["align_threads"],
        mem_gb = 32,
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
            --UNMAP_CONTAMINANT_READS true {params.daragon_workflow_params}
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
        mem_mb=lambda wildcards, attempt: (4 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 4  * attempt,
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
        mem_gb=lambda wildcards, attempt: 16  * attempt,
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
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
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

    threads: 4
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
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
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

    threads: 4
    benchmark: "benchamrks/BaseRecalibrator/{sample}.txt"
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
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

    shell:
        """
        bedtools slop -i {input.bed} -g {params.ref}.fai -b {params.padding} > {output}
        """

rule GenerateSubsettedContaminationResources:
    input: 
        intervalslist = "resources/padded.bed",
        UD = "/home/marc/Desktop/data/refs/GATK-Resources/Homo_sapiens_assembly38.contam.UD",
        BED = "/home/marc/Desktop/data/refs/GATK-Resources/Homo_sapiens_assembly38.contam.bed",
        MU = "/home/marc/Desktop/data/refs/GATK-Resources/Homo_sapiens_assembly38.contam.mu",

    conda: "../env/wes_gatk.yml"

    benchmark: "benchamrks/GenSubCont/GenerateSubsettedContaminationResources.txt"

    output: 
        UD = "03-1_SubsettedContamination/EXOME_Contams.UD",
        BED = "03-1_SubsettedContamination/EXOME_Contams.bed",
        MU = "03-1_SubsettedContamination/EXOME_Contams.mu",
        target_overlap_counts = "03-1_SubsettedContamination/target_overlap_counts.txt"

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

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

## from https://github.com/broadinstitute/warp/blob/develop/tasks/broad/BamProcessing.wdl
# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling

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

    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        '''
        verifybamid2 \
            --Verbose \
            --NumPC 4 \
            --NumThread {threads} \
            --Output {params.prefix} \
            --BamFile {input.target_bam} \
            --Reference {input.ref} \
            --SVDPrefix {params.svdprefix} \

        # extract the contamination value from the selfSM file  
        python3 {params.get_con_file} \
            --input {output.selfSM} \
            --output {output.contamination} 

        '''

#-----------------------------------------------#
#                    BAM QC                     #
#-----------------------------------------------#

# Input is the final processed bam file

# Collect sequencing yield quality metrics
rule CollecQualityYieldMetrics:
    input:
        "03_bamPrep/{sample}.bam",
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
        "03_bamPrep/{sample}.bam",
    output:
        directory("03-3_bamQC/CBQmatrix/{sample}"),
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
            --PROGRAM CollectInsertSizeMetrics \
            --PROGRAM CollectSequencingArtifactMetrics \
            --PROGRAM QualityScoreDistribution \
            --PROGRAM CollectGcBiasMetrics \
            --PROGRAM CollectBaseDistributionByCycle \
            --METRIC_ACCUMULATION_LEVEL null \
            --METRIC_ACCUMULATION_LEVEL SAMPLE \
            --METRIC_ACCUMULATION_LEVEL LIBRARY \
            --METRIC_ACCUMULATION_LEVEL READ_GROUP \
            --METRIC_ACCUMULATION_LEVEL ALL_READS

        """

# Collect ConvertSequencingArtifactToOxoG
rule ConvertSequencingArtifactToOxoG:
    input:
        bam = "03_bamPrep/{sample}.bam",
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
# rule CrossCheckFingerprints:
#     input:
#         bam = "03_bamPrep/{sample}.bam",
#         pre_adapter_detail_metrics = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix.pre_adapter_detail_metrics",
#     output:
#         "03-3_bamQC/SequencingArtifactToOxoG/{sample}.oxog_metrics"

#     conda: "../env/wes_gatk.yml"
#     params: 
#         prefix = "03-3_bamQC/SequencingArtifactToOxoG/{sample}",
#         inPasename = "03-3_bamQC/CBQmatrix/{sample}/BamMatrix",
#         ref = ref_fasta,
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt
#     threads: 1
#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             ConvertSequencingArtifactToOxoG \
#             --INPUT_BASE {params.inPasename} \
#             --OUTPUT_BASE {params.prefix} \
#             --REFERENCE_SEQUENCE {params.ref} \
#         """

#################################################
#             07 variant callign                #
#################################################

# TODO:  update once defined this rule!
Dragon_STR_Model_path = ""

rule HaplotypeCaller:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
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
        runtime = lambda wildcards, attempt: 60 * 2 * attempt,
        reduced = lambda wildcards, attempt: int(attempt * (12 * 0.80))

    threads: 4

    shell:
        """
        gatk --java-options "-Xmx{resources.reduced}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            HaplotypeCaller -R {params.ref} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation {params.extra_args} \
            -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            --interval-padding {params.padding} \
            -I {input} --native-pair-hmm-threads {threads} -ERC GVCF -O {output.vcf} \
            -contamination 0 {params.Dragon_mode} \
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

rule Annovar:
    input:
        "04_calling/{type}/variants_genotyped.filttered.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/Annovar/{type}/variants_genotyped_annotation.txt"
    output:
        directory("05_Annotation/ANNOVAR/{type}")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
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