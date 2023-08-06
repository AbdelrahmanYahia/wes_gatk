
rule HaplotypeCaller:
    input: "03_bamPrep/{sample}.pqsr.bam"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/{sample}_raw.gvcf.gz"
    params: 
        ref = ref_fasta,
        bed = bed_file,
        extra_args = config["caller_extra_args"],
        padding = config["padding"],
        WRAP_args35 = "--max_alternate_alleles 3 -variant_index_parameter 128000 -variant_index_type LINEAR -contamination 0 --read_filter OverclippedRead",
        # WRAP_args4 = "-contamination 0 ",
        # DragenMode = lambda wildcards: f"--dragen-mode" if config['use_dragon_mod'] else "",
        # dontUse_SEG = f"--disable-spanning-event-genotyping" if config['dont_use_SEG'] else "",


        # ~{if defined(dragstr_model) then "--dragstr-params-path " + dragstr_model else ""} 
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
            --floor-blocks -OVI \
            -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -L {params.bed} \
            -do-qual-approx \
            --interval-padding {params.padding} \
            -V {input} -O {output}

        """

"""
hard filter GVCF:
    gatk --java-options "-Xms2000m -Xmx2500m" \
      VariantFiltration \
      -V ~{input_vcf} \
      -L ~{interval_list} \
      --filter-expression "QD < 2.0 || FS > 30.0 || SOR > 3.0 || MQ < 40.0 || MQRankSum < -3.0 || ReadPosRankSum < -3.0" \
      --filter-name "HardFiltered" \
      -O ~{output_vcf_name}

DragenHardFilterVcf 
     gatk --java-options "-Xms2000m -Xmx2500m" \
      VariantFiltration \
      -V ~{input_vcf} \
      --filter-expression "QUAL < 10.4139" \
      --filter-name "DRAGENHardQUAL" \
      -O ~{output_vcf_name}

CNNScoreVariants 
     gatk --java-options "-Xmx10000m" CNNScoreVariants \
       -V ~{input_vcf} \
       -R ~{ref_fasta} \
       -O ~{output_vcf} \
       ~{bamout_param} \
       -tensor-type ~{tensor_type}

FilterVariantTranches 
    gatk --java-options "-Xmx6000m" FilterVariantTranches \
      -V ~{input_vcf} \
      -O ~{vcf_basename}.filtered.vcf.gz \
      ~{sep=" " prefix("--snp-tranche ", snp_tranches)} \
      ~{sep=" " prefix("--indel-tranche ", indel_tranches)} \
      --resource ~{hapmap_resource_vcf} \
      --resource ~{omni_resource_vcf} \
      --resource ~{one_thousand_genomes_resource_vcf} \
      --resource ~{dbsnp_resource_vcf} \
      --info-key ~{info_key} \
      --create-output-variant-index true

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
        # use_interval = lambda wildcards: f"-L {bed_file}" if config['use_supplied_interval'] else "",


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

"""
ScatterIntervalList 
    set -e
    mkdir out
    java -Xms1000m -Xmx1500m -jar /usr/gitc/picard.jar \
      IntervalListTools \
      SCATTER_COUNT=~{scatter_count} \
      SUBDIVISION_MODE=BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      UNIQUE=true \
      SORT=true \
      BREAK_BANDS_AT_MULTIPLES_OF=~{break_bands_at_multiples_of} \
      INPUT=~{interval_list} \
      OUTPUT=out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
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
        infile = "04-2_GenomicsDB/Interval_{interval}",
        intervals = "resources/Scatterred_intervals/00{interval}-scattered.interval_list"
    
    conda: "../env/wes_gatk.yml"
    benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

    output:
        "04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz"
    params:
        ref = ref_fasta,
        dbsnp = "~/Desktop/data/refs/vcf/resources_broad_hg38_v0_1000G_omni2.5.hg38.vcf.gz",
        Additional_annotation = ""

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
            -R {params.ref} -V gendb://{input.infile} -O {output} \
            -D {params.dbsnp} \
            --only-output-calls-starting-in-intervals \
            {params.Additional_annotation} --merge-input-intervals \
            -L {input.intervals} \

        """


## won't use as it is not the default
# rule GnarlyGenotyper:
#     input:
#         in = "04-2_GenomicsDB/Interval_{interval}",
#         intervals = "resources/Processed_intervals.bed"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         "04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz"
#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = ""
#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             GnarlyGenotyper \
#             --only-output-calls-starting-in-intervals \
#             -D {params.dbsnp} \
#             -L {input.intervals} \
#             -stand-call-conf 10 \
#             --max-alternate-alleles 5 \
#             --merge-input-intervals \
#             -R {params.ref} -V gendb://{input.in} -O {output}
#         """


# rule HardFilterAndMakeSitesOnlyVcf:
#     input:
#         "04-3_Genotyping/Scatterred/Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             VariantFiltration \
#             --filter-expression "ExcessHet > 54.69" \
#             --filter-name ExcessHet \
#             -O {output.hardfiltered} \
#             -V {input}

#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             MakeSitesOnlyVcf \
#             -I {input} \
#             -O {output.siteOnly}

#         """

# rule IndelsVariantRecalibrator:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             VariantRecalibrator \
#             -V ~{sites_only_variant_filtered_vcf} \
#             -O ~{recalibration_filename} \
#             --tranches-file ~{tranches_filename} \
#             --trust-all-polymorphic \
#             -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
#             -an ~{sep=' -an ' recalibration_annotation_values} \
#             ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
#             -mode INDEL \
#             --max-gaussians ~{max_gaussians} \
#             -resource:mills,known=false,training=true,truth=true,prior=12 ~{mills_resource_vcf} \
#             -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ~{axiomPoly_resource_vcf} \
#             -resource:dbsnp,known=true,training=false,truth=false,prior=2 ~{dbsnp_resource_vcf}

#         """


# rule SNPsVariantRecalibratorCreateModel:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             VariantRecalibrator \
#             -V ~{sites_only_variant_filtered_vcf} \
#             -O ~{recalibration_filename} \
#             --tranches-file ~{tranches_filename} \
#             --trust-all-polymorphic \
#             -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
#             -an ~{sep=' -an ' recalibration_annotation_values} \
#             ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
#             -mode SNP \
#             --sample-every-Nth-variant ~{downsampleFactor} \
#             --output-model ~{model_report_filename} \
#             --max-gaussians ~{max_gaussians} \
#             -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
#             -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
#             -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
#             -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}
#         """


# rule SNPsVariantRecalibrator:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#         MODEL_REPORT=~{model_report}

#         gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
#             VariantRecalibrator \
#             -V ~{sites_only_variant_filtered_vcf} \
#             -O ~{recalibration_filename} \
#             --tranches-file ~{tranches_filename} \
#             --trust-all-polymorphic \
#             -tranche ~{sep=' -tranche ' recalibration_tranche_values} \
#             -an ~{sep=' -an ' recalibration_annotation_values} \
#             ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
#             -mode SNP \
#             ~{model_report_arg} \
#             --max-gaussians ~{max_gaussians} \
#             -resource:hapmap,known=false,training=true,truth=true,prior=15 ~{hapmap_resource_vcf} \
#             -resource:omni,known=false,training=true,truth=true,prior=12 ~{omni_resource_vcf} \
#             -resource:1000G,known=false,training=true,truth=false,prior=10 ~{one_thousand_genomes_resource_vcf} \
#             -resource:dbsnp,known=true,training=false,truth=false,prior=7 ~{dbsnp_resource_vcf}

#         """


# rule GatherTranches:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#     set -euo pipefail

#     tranches_fofn=~{write_lines(tranches)}

#     # Jose says:
#     # Cromwell will fall over if we have it try to localize tens of thousands of files,
#     # so we manually localize files using gsutil.
#     # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
#     # PAPI doesn't do.

#     # This is here to deal with the JES bug where commands may be run twice
#     rm -rf tranches
#     mkdir tranches
#     RETRY_LIMIT=5

#     count=0
#     until cat $tranches_fofn | gsutil -m cp -L cp.log -c -I tranches/; do
#       sleep 1
#       ((count++)) && ((count >= $RETRY_LIMIT)) && break
#     done
#     if [ "$count" -ge "$RETRY_LIMIT" ]; then
#       echo 'Could not copy all the tranches from the cloud' && exit 1
#     fi

#     cat $tranches_fofn | rev | cut -d '/' -f 1 | rev | awk '{print "tranches/" $1}' > inputs.list

#     gatk --java-options "-Xmx6000m -Xmx7000m" \
#       GatherTranches \
#       --input inputs.list \
#       --mode ~{mode} \
#       --output ~{output_filename}
#         """

# rule ApplyRecalibration:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     gatk --java-options "-Xms5000m -Xmx6500m" \
#       ApplyVQSR \
#       -O tmp.indel.recalibrated.vcf \
#       -V ~{input_vcf} \
#       --recal-file ~{indels_recalibration} \
#       ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
#       --tranches-file ~{indels_tranches} \
#       --truth-sensitivity-filter-level ~{indel_filter_level} \
#       --create-output-variant-index true \
#       -mode INDEL

#     gatk --java-options "-Xms5000m -Xmx6500m" \
#       ApplyVQSR \
#       -O ~{recalibrated_vcf_filename} \
#       -V tmp.indel.recalibrated.vcf \
#       --recal-file ~{snps_recalibration} \
#       ~{true='--use-allele-specific-annotations' false='' use_allele_specific_annotations} \
#       --tranches-file ~{snps_tranches} \
#       --truth-sensitivity-filter-level ~{snp_filter_level} \
#       --create-output-variant-index true \
#       -mode SNP

#         """


# rule GatherVcfs:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     set -euo pipefail

#     # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
#     # This argument disables expensive checks that the file headers contain the same set of
#     # genotyped samples and that files are in order by position of first record.
#     gatk --java-options "-Xms6000m -Xmx6500m" \
#       GatherVcfsCloud \
#       --ignore-safety-checks \
#       --gather-type BLOCK \
#       --input ~{sep=" --input " input_vcfs} \
#       --output ~{output_vcf_name}

#     tabix ~{output_vcf_name}

      
#         """



# rule SelectFingerprintSiteVariants:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     set -euo pipefail

#     function hdb_to_interval_list() {
#       input=$1
#       awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
#     }

#     hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

#     gatk --java-options "-Xms6000m -Xmx7000m" \
#       SelectVariants \
#       --variant ~{input_vcf} \
#       --intervals hdb.interval_list \
#       --output ~{base_output_name}.vcf.gz
      
#         """


# rule CollectVariantCallingMetrics:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     gatk --java-options "-Xms6000m -Xmx7000m" \
#       CollectVariantCallingMetrics \
#       --INPUT ~{input_vcf} \
#       --DBSNP ~{dbsnp_vcf} \
#       --SEQUENCE_DICTIONARY ~{ref_dict} \
#       --OUTPUT ~{metrics_filename_prefix} \
#       --THREAD_COUNT 8 \
#       --TARGET_INTERVALS ~{interval_list}
      
#         """


# rule GatherVariantCallingMetrics:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     input_details_fofn=~{write_lines(input_details)}
#     input_summaries_fofn=~{write_lines(input_summaries)}

#     # Jose says:
#     # Cromwell will fall over if we have it try to localize tens of thousands of files,
#     # so we manually localize files using gsutil.
#     # Using gsutil also lets us parallelize the localization, which (as far as we can tell)
#     # PAPI doesn't do.

#     # This is here to deal with the JES bug where commands may be run twice
#     rm -rf metrics

#     mkdir metrics
#     RETRY_LIMIT=5

#     count=0
#     until cat $input_details_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
#       sleep 1
#       ((count++)) && ((count >= $RETRY_LIMIT)) && break
#     done
#     if [ "$count" -ge "$RETRY_LIMIT" ]; then
#       echo 'Could not copy all the metrics from the cloud' && exit 1
#     fi

#     count=0
#     until cat $input_summaries_fofn | gsutil -m cp -L cp.log -c -I metrics/; do
#       sleep 1
#       ((count++)) && ((count >= $RETRY_LIMIT)) && break
#     done
#     if [ "$count" -ge "$RETRY_LIMIT" ]; then
#       echo 'Could not copy all the metrics from the cloud' && exit 1
#     fi

#     INPUT=$(cat $input_details_fofn | rev | cut -d '/' -f 1 | rev | sed s/.variant_calling_detail_metrics//g | awk '{printf("--INPUT metrics/%s ", $1)}')

#     gatk --java-options "-Xms2000m -Xmx2500m" \
#       AccumulateVariantCallingMetrics \
#       $INPUT \
#       --OUTPUT ~{output_prefix}
#         """


# rule CrossCheckFingerprint:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#      gvcfInputsList=~{write_lines(gvcf_paths)}
#     vcfInputsList=~{write_lines(vcf_paths)}

#     cp $gvcfInputsList gvcf_inputs.list
#     cp $vcfInputsList vcf_inputs.list

#     gatk --java-options "-Xms~{java_mem}m -Xmx~{java_mem}m" \
#       CrosscheckFingerprints \
#       --INPUT gvcf_inputs.list \
#       --SECOND_INPUT vcf_inputs.list \
#       --HAPLOTYPE_MAP ~{haplotype_database} \
#       --INPUT_SAMPLE_FILE_MAP ~{sample_name_map} \
#       --CROSSCHECK_BY SAMPLE \
#       --CROSSCHECK_MODE CHECK_SAME_SAMPLE \
#       --NUM_THREADS ~{cpu} \
#       ~{true='--EXIT_CODE_WHEN_MISMATCH 0' false='' scattered} \
#       --OUTPUT ~{output_name}

#     if ~{scattered}; then
#       # UNEXPECTED_MATCH is not possible with CHECK_SAME_SAMPLE
#       matches=$(grep "EXPECTED_MATCH" ~{output_name} | wc -l)

#       # check inconclusive samples
#       expectedInconclusiveSamples=("~{sep='" "' expected_inconclusive_samples}")
#       inconclusiveSamplesCount=0
#       inconclusiveSamples=($(grep 'INCONCLUSIVE' ~{output_name} | cut -f 1))
#       for sample in ${inconclusiveSamples[@]}; do
#       if printf '%s\n' ${expectedInconclusiveSamples[@]} | grep -P '^'${sample}'$'; then
#         inconclusiveSamplesCount=$((inconclusiveSamplesCount+1))
#       fi
#       done

#       total_matches=$((inconclusiveSamplesCount + matches))
#       if [[ ${total_matches} -eq ~{num_gvcfs} ]]; then
#         >&2 echo "Found the correct number of matches (~{num_gvcfs}) for this shard"
#       else
#         >&2 echo "ERROR: Found $total_matches 'EXPECTED_MATCH' records, but expected ~{num_gvcfs}"
#       exit 1
#       fi
#     fi
#         """



# rule GatherPicardMetrics:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     # Don't use this task to gather tens of thousands of files.
#     # Cromwell can't handle it.

#     # This cannot gather metrics with histograms

#     head -n 7 ~{metrics_files[0]} > ~{output_file_name}

#     for metrics_file in ~{sep=' ' metrics_files}; do
#       sed -n '1,7d;p' $metrics_file | grep -v '^$' >> ~{output_file_name}
#     done

#         """



# rule GetFingerprintingIntervalIndices:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     function rename_intervals(){
#       interval_list=$1
#       name=$2

#       awk 'BEGIN{FS=IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {$5="'$name'"; print}' $interval_list
#     }
#     export -f rename_intervals

#     function hdb_to_interval_list(){
#       input=$1

#       awk 'BEGIN{IFS="\t";OFS="\t";} $0~"^@"{print;next;} $0~"#CHROM"{next;} {print $1,$2,$2,"+","interval-"NR}' $1
#     }

#     function rename_scatter(){
#       file=$1
#       number=$(echo $file | sed -E 's|([0-9]+)-scattered\.interval.*|\1|')
#       rename_intervals $file $number > scattered.renamed.$number.interval_list
#     }
#     export -f rename_scatter

#     # rename the intervals within each interval_list according to the number in the name of the list

#     cp ~{sep=' ' unpadded_intervals} ./

#     cat ~{write_lines(unpadded_intervals)} | xargs -n1 basename | xargs -I{} bash -c 'rename_scatter $@' _ {}

#     #find the first header
#     find . -name "scattered.renamed.*.interval_list" | head -n1 | xargs cat | grep '^@' > all.interval_list

#     # concatenate the resulting intervals (with no header)
#     find . -name "scattered.renamed.*.interval_list"  | xargs cat | grep -v '^@' >> all.interval_list

#     # convert the Haplotype_database to an interval_list
#     hdb_to_interval_list ~{haplotype_database} > hdb.interval_list

#     # find the intervals that overlap the haplotype_database
#     gatk --java-options "-Xms3000m -Xmx3250m" \
#     IntervalListTools \
#       -ACTION OVERLAPS \
#       -O all.sorted.interval_list \
#       -I all.interval_list \
#       -SI hdb.interval_list

#     if grep -v '^@' all.sorted.interval_list; then
#       grep -v '^@' all.sorted.interval_list | awk '{FS="\t"; print $5}' | uniq > indices.out
#     else
#       touch indices.out
#     fi
#         """



# rule PartitionSampleNameMap:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """

#     cut -f 2 ~{sample_name_map} > sample_paths
#     split -l ~{line_limit} -d sample_paths partition_

#     # Let the OS catch up with creation of files for glob command
#     sleep 1

#         """


# rule CalculateAverageAnnotations:
#     input:
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"
    
#     conda: "../env/wes_gatk.yml"
#     benchmark: "benchamrks/GenotypeGvcfs/Interval_{interval}.txt"

#     output:
#         hardfiltered = "04-3_Genotyping/01_Filterration/HardFilter_Interval_{interval}.vcf.gz",
#         siteOnly = "04-3_Genotyping/01_Filterration/siteOnly_Interval_{interval}.vcf.gz"

#     params:
#         ref = ref_fasta,
#         dbsnp = "",
#         Additional_annotation = "",
#         excess_het_threshold = 54.69

#     threads: 4
#     resources:
#         mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
#         # cores=config["general_low_threads"],
#         mem_gb=lambda wildcards, attempt: 8  * attempt,
#         # nodes = 1,
#         runtime = lambda wildcards, attempt: 60 * 2 * attempt

#     shell:
#         """
#     gatk --java-options "-Xms~{memory_mb-2000}m" \
#       CalculateAverageCombinedAnnotations \
#       -V ~{vcf} \
#       --summed-annotation-to-divide ~{sep=" --summed-annotation-to-divide " annotations_to_divide} \
#       -O ~{basename}.avg.vcf.gz

#         """


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
        "04_calling/variants_genotyped.gvcf.gz"

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
        gatk --java-options '-DGATK_STACKTRACE_ON_USER_EXCEPTION=true' \
            CombineGVCFs \
            -R {params.ref} \
            {params.gvcfs} \
            -G StandardAnnotation -G StandardHCAnnotation \
            -G AS_StandardAnnotation \
            --allow-old-rms-mapping-quality-annotation-data \
            -O {output}
        """
