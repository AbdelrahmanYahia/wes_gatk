
rule FastqToSam:
    input:
        R1 = f"{samples_dir}/{{sample}}_{{unit}}1.{EXT}",
        R2 = f"{samples_dir}/{{sample}}_{{unit}}2.{EXT}"
    
    conda: "../env/wes_gatk.yml"

    output:
        ubam = "0_samples/{sample}/{sample}_{unit}.ubam"


    benchmark: "benchamrks/QC/{sample}/{sample}_{unit}_trim.txt"
    threads: 4
    params:
        size = 4,
        quality = 10,
        extra = '',
        minlen = 50,
        ext = EXT

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

        ext={params.ext}
        R1={input.R1}
        echo "R1 $R1"
        SM={wildcards.sample}
        echo "SM $SM"
        PL="Illumina"
        LB=$SM
        echo "LB $LB"      
        if [[ $ext == *.gz ]]; then
            RGID=$(zcat {input.R1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        else
            RGID=$(head {input.R1} -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        fi
        echo "RGID $RGID"
        ## TODO: confirm to use this
        PU=$RGID.$LB 
        echo "PU $PU" 

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