# CNV GATK WORKFLOW 
# based on: https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants

# resources
## useful snakemake tutorial: https://evodify.com/gatk-cnv-snakemake/
## GATK github cwl file: https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/germline/cnv_germline_cohort_workflow.wdl

######## 1. Collect raw counts data with PreprocessIntervals and CollectReadCounts

# For exome data, pad target regions
rule PreprocessIntervals:
    # PreprocessIntervals pads exome targets and bins WGS intervals. 
    # Binning refers to creating equally sized intervals across the reference. 
    input:
        ""
    output:
        "targets.preprocessed.interval_list"
    params:
        ref = ref_fasta,
        bed = bed_file,

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # padding 250bp as in GATKs tutorial
        # --bin-length 0 to disable binning
        gatk PreprocessIntervals \
            -R {params.ref} \
            -L {params.bed} \
            --padding 250 \
            --bin-length 0 \
            -imr OVERLAPPING_ONLY \
            -O targets.preprocessed.interval_list
        """


"""
The tutorial generates text-based TSV (tab-separated-value) format data instead 
of the default HDF5 format by adding --format TSV to the command.
Omit this option to generate the default HDF5 format. Downstream tools process 
HDF5 format more efficiently.
Here and elsewhere in the workflow, set --interval-merging-rule (-imr) to OVERLAPPING_ONLY,
to prevent the tool from merging abutting intervals.
The tool employs a number of engine-level read filters. Of note are NotDuplicateReadFilter
and MappingQualityReadFilter. This means the tool excludes reads marked as duplicate 
and excludes reads with mapping quality less than 10. Change the mapping 
quality threshold with the --minimum-mapping-quality option.
"""
rule CollectReadCounts:
    # CollectReadCounts tabulates the raw integer counts of reads overlapping an interval.
    input:
        intervals = "targets.preprocessed.interval_list",
        sample = ""
    output:
        "sample.hdf5"
    params:
        ref = ref_fasta,
        bed = bed_file,

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # Here and elsewhere in the workflow, set --interval-merging-rule (-imr) to OVERLAPPING_ONLY, to prevent the tool from merging abutting intervals.
        # The tool employs a number of engine-level read filters. Of note are NotDuplicateReadFilter and MappingQualityReadFilter. This means the tool excludes reads marked as duplicate and excludes reads with mapping quality less than 10. Change the mapping quality threshold with the --minimum-mapping-quality option.

        # OVERLAPPING_ONLY prevents the merging of abutting intervals as recommended by the GATK team.

        gatk CollectReadCounts \
            -R {params.ref} \
            -L {input.intervals} \
            -imr OVERLAPPING_ONLY \
            -I {input.sample} \
            -O {output.sample}
        """
    

############ 2. Annotate intervals with features and subset regions of interest with FilterIntervals
"""
Explicit GC-correction, although optional, is recommended. 
The default v4.1.0.0 cnv_germline_cohort_workflow.wdl pipeline 
workflow omits explicit gc-correction and we activate it in the pipeline 
by setting do_explicit_gc_correction":"True". 
The tutorial illustrates the optional AnnotateIntervals step by performing the recommended explicit GC-content-based filtering.
"""
# You can annotate intervals with GC content, mappability, and segmental duplication information (the last two is done by adding the files as descriped below)
rule AnnotateIntervals:
    #Towards deciding which regions to exclude, AnnotateIntervals labels the given intervals with GC content and additionally with mappability and segmental duplication content if given the respective optional resource files.
    input:
        intervals = "targets.preprocessed.interval_list"
    output:
        "annotated.ref.tsv"
    params:
        ref = ref_fasta

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # The tool requires the -R reference and the -L intervals. The tool calculates GC-content for the intervals using the reference.
        # TODO: Although optional for the tool, we recommend annotating mappability by providing a --mappability-track regions file in either .bed or .bed.gz format. Be sure to merge any overlapping intervals beforehand. The tutorial omits use of this resource.
        # TODO: Optionally and additionally, annotate segmental duplication content by providing a --segmental-duplication-track regions file in either .bed or .bed.gz format

        gatk AnnotateIntervals \
            -L {input.intervals} \
            -R {params.ref} \
            -imr OVERLAPPING_ONLY \
            -O {output}
        """

# Annotated intervals are then filtered based on tunable thresholds:
# FilterIntervals takes preprocessed intervals and either annotated intervals or read counts or both.
rule FilterIntervals:
    #Towards deciding which regions to exclude, AnnotateIntervals labels the given intervals with GC content and additionally with mappability and segmental duplication content if given the respective optional resource files.
    input:
        annotatedIntervals = "annotated.ref.tsv",
        preprocessedIntervals = "targets.preprocessed.interval_list",
        all_samples = ["-I samples and repeat"]
        # TODO: double check which samples to include?
        # I guess in cohort mode we chose most or all samples! 

    output:
        " ref.cohort.gc.filtered.interval_list"
    params:
        ref = ref_fasta,
        samples = lambda:return samples

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # The tool requires the -R reference and the -L intervals. The tool calculates GC-content for the intervals using the reference.
        # TODO: Although optional for the tool, we recommend annotating mappability by providing a --mappability-track regions file in either .bed or .bed.gz format. Be sure to merge any overlapping intervals beforehand. The tutorial omits use of this resource.
        # TODO: Optionally and additionally, annotate segmental duplication content by providing a --segmental-duplication-track regions file in either .bed or .bed.gz format
        # TODO: To follow gatk github repo change:
        #       --low-count-filter-percentage-of-samples 100 
        #       --extreme-count-filter-percentage-of-samples 100 

        gatk FilterIntervals \
            -L {input.preprocessedIntervals}\
            --annotated-intervals {input.annotatedIntervals} \
            -I {params.samples} !samples-as-in-mergebams! \
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
            --interval-merging-rule OVERLAPPING_ONLY \
            -O {output}
        """


############ 3. Call autosomal and allosomal contig ploidy with DetermineGermlineContigPloidy
# calling contig ploidy is needed to generate global baseline coverage and noise data for the subsequent steps:
## DetermineGermlineContigPloidy calls contig level ploidies for both autosomal, e.g. human chr20, and allosomal contigs, e.g. human chrX. The tool determines baseline contig ploidies using sample coverages and contig ploidy priors that give the prior probabilities for each ploidy state for each contig. In this process, the tool generates global baseline coverage and noise data GermlineCNVCaller will use later
rule DetermineGermlineContigPloidy:
    #Towards deciding which regions to exclude, AnnotateIntervals labels the given intervals with GC content and additionally with mappability and segmental duplication content if given the respective optional resource files.
    input:
        filteredIntervals = "ref.cohort.gc.filtered.interval_list",
        all_samples_intervals = ["-I samples and repeat"]
        # TODO: double check which samples to include?

    output:
        " ref.cohort.gc.filtered.interval_list"
    params:
        ref = ref_fasta,
        samples = lambda:return samples,
        contigPloidyPriors = "path/to/tsv"

    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # This produces a ploidy-calls directory and a ploidy-model directory. 
        # The ploidy-calls directory contains a folder of data for each sample 
        # in the cohort including the contig ploidy calls. Each sample directory, 
        # e.g. ploidy-calls/SAMPLE_0, contains five files.

        # The ploidy-model directory contains aggregated model data for the cohort.
        # This is the model to provide to a case-mode DetermineGermlineContigPloidy 
        # analysis and to GermlineCNVCaller. The tutorial ploidy-model directory 
        # contains the eight files as follows.

        # The cohort mode requires a --contig-ploidy-priors table and produces a ploidy model.
        gatk DetermineGermlineContigPloidy  \
            -L {input.filteredIntervals}\
            -I {params.all_samples_intervals} !samples-as-in-mergebams! \
            -imr OVERLAPPING_ONLY \
            --contig-ploidy-priors {params.contigPloidyPriors} \
            -O {output.dir} \
            --output-prefix ploidy \
            --verbosity DEBUG
        """



############ 4. Call copy number variants with GermlineCNVCaller
# GermlineCNVCaller learns a denoising model per scattered 
# shard while consistently calling CNVs across the shards. 
# The tool models systematic biases and CNVs simultaneously,
#  which allows for sensitive detection of both rare and common CNVs

# GermlineCNVCaller in COHORT MODE
# Call gCNVs on the 24-sample cohort in two scatters. Notice the different -L intervals and --output-prefix basenames.


## TODO: split intervals for parallelization
## follow this: https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants#4.2

rule GermlineCNVCaller:
    input:


    output:
    params:


    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        # TODO: modify this to work with this workflow
        # This produces per-interval gCNV calls for each of the cohort samples and a gCNV model for the cohort. 

        # How do I increase the sensitivity of detection?
        ##   For tuning, first consider the coherence length 
        ##   parameters, p-alt, p-active and the psi-scale parameters. 
        ##   These hyperparameters are just a few of the plethora of adjustable 
        ##   parameters GermlineCNVCaller offers. Refer to the GermlineCNVCaller 
        ##   tool documentation for detailed explanations, and ask on the 
        ##   GATK Forum for further guidance.

        # __WGS parameters that increase the sensitivity of calling from Mark Walker
        ###     --class-coherence-length 1000.0 \
        ###     --cnv-coherence-length 1000.0 \
        ###     --enable-bias-factors false \
        ###     --interval-psi-scale 1.0E-6 \
        ###     --log-mean-bias-standard-deviation 0.01 \
        ###     --sample-psi-scale 1.0E-6 \

        gatk GermlineCNVCaller \
            --run-mode COHORT \
            -L scatter-sm/twelve_1of2.interval_list \
            -I cvg/HG00096.tsv -I cvg/HG00268.tsv -I cvg/HG00419.tsv -I cvg/HG00759.tsv \
            -I cvg/HG01051.tsv -I cvg/HG01112.tsv -I cvg/HG01500.tsv -I cvg/HG01565.tsv \
            -I cvg/HG01583.tsv -I cvg/HG01595.tsv -I cvg/HG01879.tsv -I cvg/HG02568.tsv \
            -I cvg/HG02922.tsv -I cvg/HG03006.tsv -I cvg/HG03052.tsv -I cvg/HG03642.tsv \
            -I cvg/HG03742.tsv -I cvg/NA18525.tsv -I cvg/NA18939.tsv -I cvg/NA19017.tsv \
            -I cvg/NA19625.tsv -I cvg/NA19648.tsv -I cvg/NA20502.tsv -I cvg/NA20845.tsv \
            --contig-ploidy-calls ploidy-calls \
            --annotated-intervals twelveregions.annotated.tsv \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output cohort24-twelve \
            --output-prefix cohort24-twelve_1of2 \
            --verbosity DEBUG
        """


"""
Comments on select sensitivity parameters

Decreasing --class-coherence-length from its default of 10,000bp to 1000bp decreases the expected length of contiguous segments. Factor for bin size when tuning.
Decreasing --cnv-coherence-length from its default 10,000bp to 1000bp decreases the expected length of CNV events. Factor for bin size when tuning.
Turning off --enable-bias-factors from the default true state to false turns off active discovery of learnable bias factors. This should always be on for targeted exome data and in general can be turned off for WGS data.
Decreasing --interval-psi-scale from its default of 0.001 to 1.0E-6 reduces the scale the tool considers normal in per-interval noise.
Decreasing --log-mean-bias-standard-deviation from its default of 0.1 to 0.01 reduces what is considered normal noise in bias factors.
Decreasing --sample-psi-scale from its default of 0.0001 to 1.0E-6 reduces the scale that is considered normal in sample-to-sample variance.


Additional parameters to consider include --depth-correction-tau, --p-active and --p-alt.

--depth-correction-tau has a default of 10000.0 (10K) and defines the precision of read-depth concordance with the global depth value.
--p-active has a default of 1e-2 (0.01) and defines the prior probability of common CNV states.
p-alt has a default of 1e-6 (0.000001) and defines the expected probability of CNV events (in rare CNV states).

"""


####### 5. Call copy number segments and consolidate sample results with PostprocessGermlineCNVCalls
# PostprocessGermlineCNVCalls consolidates the scattered GermlineCNVCaller results, performs segmentation and calls copy number states. The tool generates per-interval and per-segment sample calls in VCF format and runs on a single sample at a time.

## Each command generates two VCFs with indices. The genotyped-intervals VCF contains variant records for each analysis bin and therefore data covers only the interval regions. For the tutorial's small data, this gives 1400 records. The genotyped-segments VCF contains records for each contiguous copy number state segment. For the tutorial's small data, this is 30 and 31 records for cohort and case mode analyses, respectively.
## The two modes--cohort and case--give highly concordant but slightly different results for sample NA19017. The factor that explains the difference is the contribution of the sample itself to the model.

rule PostprocessGermlineCNVCalls:
    input:


    output:
    params:


    conda: "../env/wes_gatk.yml"

    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt

    shell:
        """
        gatk PostprocessGermlineCNVCalls \
            --model-shard-path cohort24-twelve/cohort24-twelve_1of2-model \
            --model-shard-path cohort24-twelve/cohort24-twelve_2of2-model \
            --calls-shard-path cohort24-twelve/cohort24-twelve_1of2-calls \
            --calls-shard-path cohort24-twelve/cohort24-twelve_2of2-calls \
            --allosomal-contig chrX --allosomal-contig chrY \
            --contig-ploidy-calls ploidy-calls \
            --sample-index 19 \
            --output-genotyped-intervals genotyped-intervals-cohort24-twelve-NA19017.vcf.gz \
            --output-genotyped-segments genotyped-segments-cohort24-twelve-NA19017.vcf.gz \
            --sequence-dictionary ref/Homo_sapiens_assembly38.dict

        """


"""
Comments on select parameters

Specify a --model-shard-path directory for each scatter of the cohort model.
Specify a --calls-shard-path directory for each scatter of the cohort or case analysis.
Specify the --contig-ploidy-calls directory for the cohort or case analysis.
By default --autosomal-ref-copy-number is set to 2.
Define allosomal contigs with the --allosomal-contig parameter.
The tool requires specifying the --output-genotyped-intervals VCF.
Optionally generate segmented VCF results with --output-genotyped-segments. The tool segments the regions between the starting bin and the ending bin on a contig. The implication of this is that even if there is a gap between two analysis bins on the same contig, if the copy number state is equal for the bins, then the bins and the entire region between can end up a part of the same segment. The extent of this smoothing over gaps depends on the --cnv-coherence-length parameter.
The --sample-index refers to the index number given to a sample by GermlineCNVCaller. In a case mode analysis of a single sample, the index will always be zero.
The --sequence-dictionary is optional. Without it, the tool generates unindexed VCF results. Alternatively, to produce the VCF indices, provide the -R reference FASTA or use IndexFeatureFile afterward. The v4.1.0.0 cnv_germline_cohort_workflow.wdl pipeline workflow omits index generation.
"""


#### TODO: implement this :
## https://gatk.broadinstitute.org/hc/en-us/articles/360035531452
## for after gCNV calling considerations