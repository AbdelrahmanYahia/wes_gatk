
rule FastqToSam:
    input:
        R1 = f"{samples_dir}/{{sample}}_{{unit}}_1.fq.gz",
        R2 = f"{samples_dir}/{{sample}}_{{unit}}_2.fq.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        ubam = "0_samples/{sample}/{sample}_{unit}.ubam"


    benchmark: "benchamrks/QC/{sample}/{sample}_{unit}_trim.txt"
    threads: 4
    params:
        size = 4,
        quality = 10,
        extra = '',
        minlen = 50

    resources:
        mem_mb=5120,
        cores=4,
        mem_gb=5,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    log: 
        "logs/trimmomatic/{sample}/{sample}_{unit}.txt"
    shell:
        '''
        R1={input.R1}
        SM={wildcards.sample}
        PL="Illumina"
        LB="{wildcards.sample}_{wildcards.unit}"
        name=$(basename $R1 | cut -d'_' -f1)
        RGID=$(head -n1 $R1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
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
            -PU $PU 
        '''

rule MarkIlluminaAdapters:
    input:
        "0_samples/{sample}/{sample}_{unit}.ubam"
    output:
        bam="0_samples/{sample}/{sample}_{unit}.adab.ubam",
        metrics="0_samples/{sample}/{sample}_{unit}.adap_metrics.txt"
    threads:4
    resources:
        mem_mb=5120,
        cores=4,
        mem_gb=5,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            MarkIlluminaAdapters \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} 
        '''