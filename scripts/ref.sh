##### https://github.com/gatk-workflows
##### ref for dev the gatk pest paractice snakemake workflow

bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
"""
-K 100000000: This parameter sets the maximum 
size of the prefix seed array in BWA. The value 
100000000 indicates that BWA will allocate memory 
for storing 100 million seed positions. Seeds are 
short sequences used to quickly identify potential 
alignment locations.

-p: This parameter is used to enable the pairing of 
reads. In DNA sequencing, paired-end reads are generated 
from both ends of a DNA fragment. The -p option tells BWA 
to consider the paired-end information during the alignment 
process.

-v 3: This parameter sets the verbosity level for BWA mem. 
The value 3 indicates a high level of verbosity, which means 
BWA will produce more detailed output during the alignment process. 
This can be useful for debugging or gaining more insights into the 
alignment procedure.
"""

### SamToFastqAndBwaMem  
# set the bash variable needed for the command-line
bash_ref_fasta=~{ref_fasta}
# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment

java -Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G -jar ~{gotc_path}picard.jar \
    SamToFastq \
    INPUT=~{input_bam} \
    FASTQ=/dev/stdout \
    INTERLEAVE=true \
    NON_PF=true \
    | \
    ~{bwa_path}~{bwa_commandline} /dev/stdin -  2> >(tee ~{output_bam_basename}.bwa.stderr.log >&2) \
    | \
    samtools view -1 - > ~{output_bam_basename}.bam

### output_bam = "~{output_bam_basename}.bam"

#################################################

### MergeBamAlignment 
# Merge original input uBAM file with BWA-aligned BAM file

    # set the bash variable needed for the command-line
    bash_ref_fasta=~{ref_fasta}
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
      MergeBamAlignment \
      --VALIDATION_STRINGENCY SILENT \
      --EXPECTED_ORIENTATIONS FR \
      --ATTRIBUTES_TO_RETAIN X0 \
      --ALIGNED_BAM ~{aligned_bam} \
      --UNMAPPED_BAM ~{unmapped_bam} \
      --OUTPUT ~{output_bam_basename}.bam \
      --REFERENCE_SEQUENCE ~{ref_fasta} \
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
      --PROGRAM_GROUP_VERSION "~{bwa_version}" \
      --PROGRAM_GROUP_COMMAND_LINE "~{bwa_commandline}" \
      --PROGRAM_GROUP_NAME "bwamem" \
      --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
      --ALIGNER_PROPER_PAIR_FLAGS true \
      --UNMAP_CONTAMINANT_READS true

### File output_bam = "~{output_bam_basename}.bam"

#######################

### SortAndFixTags 
# Sort BAM file by coordinate order and fix tag values for NM and UQ

    set -o pipefail

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_sort}G" \
      SortSam \
      --INPUT ~{input_bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
    | \
    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb_fix}G" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ~{ref_fasta}

##    File output_bam = "~{output_bam_basename}.bam"
##    File output_bam_index = "~{output_bam_basename}.bai"
##    File output_bam_md5 = "~{output_bam_basename}.bam.md5"

#######################

### MarkDuplicates 
# Mark duplicate reads to avoid counting non-independent observations

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
      MarkDuplicates \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --METRICS_FILE ~{metrics_filename} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true

    ### File output_bam = "~{output_bam_basename}.bam"
    ### File duplicate_metrics = "~{metrics_filename}"


########################


### CreateSequenceGroupingTSV 
# Generate sets of intervals for scatter-gathering over chromosomes

   python <<CODE
    with open("~{ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
    # some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
            temp_size = sequence_tuple[1]
    # add the unmapped sequences as a separate line to ensure that they are recalibrated as well
    with open("sequence_grouping.txt","w") as tsv_file:
      tsv_file.write(tsv_string)
      tsv_file.close()

    tsv_string += '\n' + "unmapped"

    with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
      tsv_file_with_unmapped.write(tsv_string)
      tsv_file_with_unmapped.close()
CODE 

 ###   Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
  ###  Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")

  ###############

### BaseRecalibrator 
# Generate Base Quality Score Recalibration (BQSR) model


    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      --use-original-qualities \
      -O ~{recalibration_report_filename} \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
      -L ~{sep=" -L " sequence_group_interval}

###     File recalibration_report = "~{recalibration_report_filename}"

##########################

### GatherBqsrReports 
# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.

    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      GatherBQSRReports \
      -I ~{sep=' -I ' input_bqsr_reports} \
      -O ~{output_report_filename}

###     File output_bqsr_report = "~{output_report_filename}"

##################

## ApplyBQSR 
# Apply Base Quality Score Recalibration (BQSR) model


    ~{gatk_path} --java-options "-Xms~{command_mem_gb}G" \
      ApplyBQSR \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{output_bam_basename}.bam \
      -L ~{sep=" -L " sequence_group_interval} \
      -bqsr ~{recalibration_report} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities

###     File recalibrated_bam = "~{output_bam_basename}.bam"

############################

### GatherBamFiles 
# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs

    ~{gatk_path} --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{command_mem_gb}G" \
      GatherBamFiles \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true

###     File output_bam = "~{output_bam_basename}.bam"
###    File output_bam_index = "~{output_bam_basename}.bai"
###    File output_bam_md5 = "~{output_bam_basename}.bam.md5"

##############################

### CramToBamTask 
    set -e
    set -o pipefail

    ~{samtools_path} view -h -T ~{ref_fasta} ~{input_cram} |
    ~{samtools_path} view -b -o ~{sample_name}.bam -
    ~{samtools_path} index -b ~{sample_name}.bam
    mv ~{sample_name}.bam.bai ~{sample_name}.bai

###    File output_bam = "~{sample_name}.bam"
###    File output_bai = "~{sample_name}.bai"

#######################

### HaplotypeCaller 
# HaplotypeCaller per-sample in GVCF mode

    set -e
  
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      HaplotypeCaller \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -L ~{interval_list} \
      -O ~{output_filename} \
      -contamination ~{default="0" contamination} \
      -G StandardAnnotation -G StandardHCAnnotation ~{true="-G AS_StandardAnnotation" false="" make_gvcf} \
      -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
      ~{true="-ERC GVCF" false="" make_gvcf} \
      ~{if defined(gcs_project_for_requester_pays) then "--gcs-project-for-requester-pays ~{gcs_project_for_requester_pays}" else ""} \
      ~{bamout_arg}

    # Cromwell doesn't like optional task outputs, so we have to touch this file.
    touch ~{vcf_basename}.bamout.bam 

###    File output_vcf = "~{output_filename}"
###    File output_vcf_index = "~{output_filename}.tbi"
###    File bamout = "~{vcf_basename}.bamout.bam"

##################

### MergeGVCFs 
# Merge GVCFs generated per-interval for the same sample


    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G"  \
      MergeVcfs \
      --INPUT ~{sep=' --INPUT ' input_vcfs} \
      --OUTPUT ~{output_filename}

###    File output_vcf = "~{output_filename}"
###    File output_vcf_index = "~{output_filename}.tbi"

#####################
##### # Basic Joint Genotyping with GATK4 (not Best Practices, just demo)

### RenameAndIndexFile 

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      IndexFeatureFile \
      -I ~{new_name} \
      -O ~{index_name}


### ImportGVCFs 

    rm -rf ~{workspace_dir_name}

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      GenomicsDBImport \
      -V ~{sep=' -V ' input_gvcfs} \
      -L ~{interval} \
      --genomicsdb-workspace-path ~{workspace_dir_name} \
      --batch-size 50 \
      --reader-threads 5 \
      --merge-input-intervals \
      --consolidate

    tar -cf ~{tarred_workspace_name} ~{workspace_dir_name}
### File output_workspace = "~{tarred_workspace_name}"
#####

### GenotypeGVCFs 
    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      GenotypeGVCFs \
      -R ~{ref_fasta} \
      -V gendb://$WORKSPACE \
      -L ~{interval} \
      -O ~{output_vcf_filename} \
      -G StandardAnnotation -G AS_StandardAnnotation \
      --allow-old-rms-mapping-quality-annotation-data \
      --merge-input-intervals

###     File output_vcf = "~{output_vcf_filename}"
###     File output_vcf_index = "~{output_vcf_filename + output_index_suffix}"


############
### MergeVCFs 

    ~{gatk_path} --java-options "-Xmx~{command_mem_gb}G ~{java_opt}" \
      MergeVcfs \
      -I ~{sep=' -I' input_vcfs} \
      -O ~{merged_vcf_filename}

###    File output_vcf = "~{merged_vcf_filename}"
###    File output_vcf_index = "~{merged_vcf_filename + output_index_suffix}"



######################### END ########################