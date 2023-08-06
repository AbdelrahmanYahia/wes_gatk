
"""
CollectQualityYieldMetrics 
    java -Xms2000m -Xmx3000m -jar /usr/picard/picard.jar \
      CollectQualityYieldMetrics \
      INPUT=~{input_bam} \
      OQ=true \
      OUTPUT=~{metrics_filename}

CollectUnsortedReadgroupBamQualityMetrics 
    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=QualityScoreDistribution \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=ALL_READS

    touch ~{output_bam_prefix}.insert_size_metrics
    touch ~{output_bam_prefix}.insert_size_histogram.pdf


CollectReadgroupBamQualityMetrics 
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics

    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=READ_GROUP

CollectAggregationMetrics 
    # These are optionally generated, but need to exist for Cromwell's sake
    touch ~{output_bam_prefix}.gc_bias.detail_metrics \
      ~{output_bam_prefix}.gc_bias.pdf \
      ~{output_bam_prefix}.gc_bias.summary_metrics \
      ~{output_bam_prefix}.insert_size_metrics \
      ~{output_bam_prefix}.insert_size_histogram.pdf

    java -Xms5000m -Xmx6500m -jar /usr/picard/picard.jar \
      CollectMultipleMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      OUTPUT=~{output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=QualityScoreDistribution \
      ~{true='PROGRAM="CollectGcBiasMetrics"' false="" collect_gc_bias_metrics} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY


ConvertSequencingArtifactToOxoG 
    input_base=$(dirname ~{pre_adapter_detail_metrics})/~{base_name}
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m \
      -jar /usr/picard/picard.jar \
      ConvertSequencingArtifactToOxoG \
      --INPUT_BASE $input_base \
      --OUTPUT_BASE ~{base_name} \
      --REFERENCE_SEQUENCE ~{ref_fasta}

CrossCheckFingerprints 
    java -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms3000m -Xmx3000m \
      -jar /usr/picard/picard.jar \
      CrosscheckFingerprints \
      OUTPUT=~{metrics_filename} \
      HAPLOTYPE_MAP=~{haplotype_database_file} \
      EXPECT_ALL_GROUPS_TO_MATCH=true \
      INPUT=~{sep=' INPUT=' input_bams} \
      LOD_THRESHOLD=~{lod_threshold} \
      CROSSCHECK_BY=~{cross_check_by}

CheckFingerprintTask 
    set -e
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -Dpicard.useLegacyParser=false -jar /usr/picard/picard.jar \
    CheckFingerprint \
      --INPUT ~{input_file} \
      ~{if defined(input_vcf) then "--OBSERVED_SAMPLE_ALIAS \"" + input_sample_alias + "\"" else ""} \
      --GENOTYPES ~{genotypes} \
      --EXPECTED_SAMPLE_ALIAS "~{expected_sample_alias}" \
      ~{if defined(input_bam) then "--IGNORE_READ_GROUPS true" else ""} \
      --HAPLOTYPE_MAP ~{haplotype_database_file} \
      --GENOTYPE_LOD_THRESHOLD ~{genotype_lod_threshold} \
      --SUMMARY_OUTPUT ~{summary_metrics_location} \
      --DETAIL_OUTPUT ~{detail_metrics_location} \
      ~{"--REFERENCE_SEQUENCE " + ref_fasta} \
      ~{true='--EXIT_CODE_WHEN_NO_VALID_CHECKS 0' false='' allow_lod_zero}

    CONTENT_LINE=$(cat ~{summary_metrics_location} |
    grep -n "## METRICS CLASS\tpicard.analysis.FingerprintingSummaryMetrics" |
    cut -f1 -d:)
    CONTENT_LINE=$(($CONTENT_LINE+2))
    sed '8q;d' ~{summary_metrics_location} | cut -f5 > lod


CheckPreValidation 
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ~{duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ~{chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python3 <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE


ValidateSamFile 
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      ValidateSamFile \
      INPUT=~{input_bam} \
      OUTPUT=~{report_filename} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      ~{"MAX_OUTPUT=" + max_output} \
      IGNORE=~{default="null" sep=" IGNORE=" ignore} \
      MODE=VERBOSE \
      ~{default='SKIP_MATE_VALIDATION=false' true='SKIP_MATE_VALIDATION=true' false='SKIP_MATE_VALIDATION=false' is_outlier_data} \
      IS_BISULFITE_SEQUENCED=false

CollectWgsMetrics 
    java -Xms2000m -Xmx2500m -jar /usr/picard/picard.jar \
      CollectWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}

CollectRawWgsMetrics 
    java -Xms~{java_memory_size}m -jar /usr/picard/picard.jar \
      CollectRawWgsMetrics \
      INPUT=~{input_bam} \
      VALIDATION_STRINGENCY=SILENT \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      INCLUDE_BQ_HISTOGRAM=true \
      INTERVALS=~{wgs_coverage_interval_list} \
      OUTPUT=~{metrics_filename} \
      USE_FAST_ALGORITHM=true \
      READ_LENGTH=~{read_length}

CollectHsMetrics 
    java -Xms~{java_memory_size}m -Xmx~{max_heap}m -jar /usr/picard/picard.jar \
      CollectHsMetrics \
      INPUT=~{input_bam} \
      REFERENCE_SEQUENCE=~{ref_fasta} \
      VALIDATION_STRINGENCY=SILENT \
      TARGET_INTERVALS=~{target_interval_list} \
      BAIT_INTERVALS=~{bait_interval_list} \
      METRIC_ACCUMULATION_LEVEL=null \
      METRIC_ACCUMULATION_LEVEL=SAMPLE \
      METRIC_ACCUMULATION_LEVEL=LIBRARY \
      OUTPUT=~{metrics_filename}

CalculateReadGroupChecksum 
    java -Xms1000m -Xmx3500m -jar /usr/picard/picard.jar \
      CalculateReadGroupChecksum \
      INPUT=~{input_bam} \
      OUTPUT=~{read_group_md5_filename}


"""