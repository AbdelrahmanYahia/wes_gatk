
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
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
#    log: 
#        "logs/fastqtosam/{sample}/{sample}-{unit}.txt"
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
        if [[ "$ext" == *".gz" ]]; then
            RGID=$(head -n1 <(zcat {input.R1}) | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
            #RGID=$(zcat {input.R1} | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
        else
            echo "non gz file"
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
            -PU $PU \
            # > {log} 2>&1
        '''

rule MarkIlluminaAdapters:
    input:
        "0_samples/{sample}/{sample}-{unit}.ubam"
    output:
        bam="0_samples/{sample}/{sample}-{unit}.adab.ubam",
        metrics="0_samples/{sample}/{sample}-{unit}.adap_metrics.txt"
    threads:4
    resources:
        mem_mb=lambda wildcards, attempt: (8 * 1024) * attempt,
        # cores=config["general_low_threads"],
        mem_gb=lambda wildcards, attempt: 8  * attempt,
        # nodes = 1,
        runtime = lambda wildcards, attempt: 60 * 2 * attempt
#   log: 
#        "logs/markilluminaAdabs/{sample}/{sample}-{unit}.txt"
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
