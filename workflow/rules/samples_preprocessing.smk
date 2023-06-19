
rule FastqToSam:
    input:
        R1 = f"{samples_dir}/{{sample}}-{{unit}}1.{EXT}",
        R2 = f"{samples_dir}/{{sample}}-{{unit}}2.{EXT}"
    
    conda: "../env/wes_gatk.yml"

    output:
        ubam = "0_samples/{sample}/{sample}-{unit}.ubam"

    threads: config["gen_ubam_threads"]
    params:
        ext = EXT

    resources:
        mem_mb = int(config["gen_ubam_mem"])*1024,
        cores = config["gen_ubam_threads"],
        mem_gb = int(config["gen_ubam_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    log: 
        "logs/fastqtosam/{sample}/{sample}-{unit}.txt"
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
        mem_mb = int(config["gen_ubam_mem"])*1024,
        cores = config["gen_ubam_threads"],
        mem_gb = int(config["gen_ubam_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    log: 
        "logs/markilluminaAdabs/{sample}/{sample}-{unit}.txt"

    shell:
        '''
        gatk --java-options "-Xmx{resources.mem_gb}G -XX:+UseParallelGC -XX:ParallelGCThreads={threads}" \
            MarkIlluminaAdapters \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            # > {log} 2>&1
        '''