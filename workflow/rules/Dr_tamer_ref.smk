rule FastqToSam:
    input:
        r1="data/fastq/{fragment}_R1_001.fastq.gz",
        r2="data/fastq/{fragment}_R2_001.fastq.gz"
    output:
        "data/uBAM/{fragment}.bam"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        name=$(basename {input.r1})
        SM=$(echo $name | cut -d "_" -f1)
        LB=$(echo $name | cut -d"_" -f1,2)  ## We use <Index.Sequence> in the Illumina file name as an index to the library
        #batch=$(basename "$(dirname {input.r1})")
        #if [ "$batch" != "trimmed" ];then LB=$batch.$LB;fi
        PL="Illumina"
        ##read Fastq 1st read, check the format.If typical, identify ID as "<instrument>:<run number>:<flowcell ID>:<lane>"
        header=$(head -n1 <(zcat {input.r1}) | grep ':*:*:*:*:*:*')
        if [ "$header" != "" ]; then
            RGID=$(echo "$header" | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        else # "make unique ID and PU using checksum"
            checksum=$(shasum {input.r1} | awk '{{ print $1 }}')
            RGID="UnChrPU_"$checksum
        fi
        PU=$RGID.$LB

        module load Java/jdk1.8.0
        source activate gatk
        #java -jar picard.jar FastqToSam
        gatk --java-options "-Xmx6G" FastqToSam \
        -F1={input.r1} \
        -F2={input.r2} \
        -O={output} \
        -SM=$SM \
        -LB=$LB \
        -PL=$PL \
        -RG=$RGID \
        -PU=$PU \
        --TMP_DIR="tmp/{wildcards.fragment}"
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.1.2/picard_illumina_MarkIlluminaAdapters.php
rule MarkIlluminaAdapters:
    input:
        "data/uBAM/{fragment}.bam"
    output:
        bam="data/adap_uBAM/{fragment}.adap.bam",
        metrics="data/adap_uBAM/{fragment}.adap_metrics.txt"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    threads:1
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        gatk --java-options "-Xmx6G" MarkIlluminaAdapters \
        -I={input} \
        -O={output.bam} \
        -M={output.metrics} \
        --TMP_DIR="tmp2/{wildcards.fragment}"
        '''

rule download_ref:
    output:
        "refGenome/canFam3_chr.fa"
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 30,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000
    shell:
        '''
        mkdir refGenome && cd refGenome
        wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/canFam3/bigZips/canFam3.fa.gz' -O canFam3.fa.gz
        gunzip canFam3.fa.gz
        cat canFam3.fa | awk '{if($1 ~ ">chrUn_"){f=0;}else if($1 ~ ">chr"){print $0;f=1;}else if(f){print $0;}}' > canFam3_chr.fa
        '''

rule bwa_index:
    input:
        "refGenome/canFam3_chr.fa"
    output:
        "refGenome/BwaIndex/genome.fa",
        expand("refGenome/BwaIndex/genome.{ext}", ext=['amb', 'ann', 'bwt', 'pac', 'sa']),
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 2,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 16
    shell:
        '''
        if [ ! -f refGenome/BwaIndex/genome.fa ];then ln -s ../canFam3_chr.fa refGenome/BwaIndex/genome.fa;fi
        module load bwa/0.7.7.r441
        bwa index -p refGenome/BwaIndex/genome -a bwtsw {input}
        '''

rule GATK_index:
    input:
        "refGenome/canFam3_chr.fa"
    output:
        ref="refGenome/gatkIndex/genome.fa",
        index="refGenome/gatkIndex/genome.fa.fai",
        dict="refGenome/gatkIndex/genome.dict",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 1,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 8
    shell:
        '''
        if [ ! -f refGenome/gatkIndex/genome.fa ];then ln -s ../canFam3_chr.fa refGenome/gatkIndex/genome.fa;fi
        module load SAMTools/1.5
        module load picardTools/1.89
        samtools faidx "refGenome/gatkIndex/genome.fa"
        java -Xmx4g -jar $PICARD/CreateSequenceDictionary.jar R= {input} O= {output.dict}
        '''

# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.2.0/picard_sam_SamToFastq.php
# https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_sam_MergeBamAlignment.php
# With the original MarkDuplicates aproach, MergeBamAlignment runs with SORT_ORDER 'coordinate'. If you run with SORT_ORDER other than 'coordinate' (e.g. in the WDL pipeline implementation), you will need to use SortSam to sort by queryname and SetNmAndUqTags (DEPRECATED: Use SetNmMdAndUqTags instead) to fix these tags
rule align:
    input:
        bam="data/adap_uBAM/{fragment}.adap.bam",
        bwa_ref="refGenome/BwaIndex/genome.fa",
        gatk_ref="refGenome/gatkIndex/genome.fa",
    output:
        bam="data/mapped_reads/{fragment}.bam",
    resources:
        walltime = lambda wildcards, attempt: 2**(attempt - 1) * 60 * 60 * 8,
        mem = lambda wildcards, attempt: 2**(attempt - 1) * 1000000000 * 32
    threads:4
    shell:
        '''
        module load Java/jdk1.8.0
        source activate gatk
        module load bwa/0.7.7.r441
        gatk --java-options "-Xmx30G" SamToFastq \
        -I={input.bam} \
        --FASTQ=/dev/stdout \
        --CLIPPING_ATTRIBUTE=XT --CLIPPING_ACTION=2 --INTERLEAVE=true -NON_PF=true \
        --TMP_DIR="tmp3/{wildcards.fragment}" | \
        bwa mem -M -t {threads} -p refGenome/BwaIndex/genome /dev/stdin | \
        gatk --java-options "-Xmx30G" MergeBamAlignment \
        --ALIGNED_BAM=/dev/stdin \
        --UNMAPPED_BAM={input.bam} \
        --OUTPUT={output.bam} \
        -R={input.gatk_ref} --CREATE_INDEX=true --ADD_MATE_CIGAR=true \
        --CLIP_ADAPTERS=false --CLIP_OVERLAPPING_READS=true \
        --INCLUDE_SECONDARY_ALIGNMENTS=true --MAX_INSERTIONS_OR_DELETIONS=-1 \
        --PRIMARY_ALIGNMENT_STRATEGY=MostDistant --ATTRIBUTES_TO_RETAIN=XS \
        --TMP_DIR="tmp4/{wildcards.fragment}"
        '''
