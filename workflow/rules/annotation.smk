
rule consequence:
    input: "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/QC/{type}/bcftools_csq.vcf"
    threads: config["general_low_threads"]
    resources:
        mem_mb=int(config["general_low_mem"])* 1024,
        cores=config["general_low_threads"],
        mem_gb=int(config["general_low_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
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

    output: "04_calling/QC/{type}/bcftools.stats"
    threads: config["general_low_threads"]
    resources:
        mem_mb=int(config["general_low_mem"])* 1024,
        cores=config["general_low_threads"],
        mem_gb=int(config["general_low_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
            bcftools stats {input} > {output}
        """

rule plot_bcftools_stats:
    input: "04_calling/QC/{type}/bcftools.stats"
    
    conda: "../env/wes_gatk.yml"

    output: directory("04_calling/QC/{type}/bcftools_plots")
    threads: config["general_low_threads"]
    resources:
        mem_mb=int(config["general_low_mem"])* 1024,
        cores=config["general_low_threads"],
        mem_gb=int(config["general_low_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    shell:
        """
        plot-vcfstats -p {output} {input}
        """

rule Nirvana:
    input:
        "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        "05_Annotation/Nirvana/{type}/Annotation.json.gz"
    threads: config["annotation_threads"]
    resources:
        mem_mb=int(config["annotation_mem"])* 1024,
        cores=config["annotation_threads"],
        mem_gb=int(config["annotation_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/",
        Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat",
        Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both",
        Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll",
        file_name = lambda wildcards: f"05_Annotation/Nirvana/{wildcards.type}/Annotation",
        extra = config["nirvana_extra_args"]

    shell:
        """
        dotnet {params.Nirvana_cmd} \
            -i {input} \
            -o {params.file_name} \
            -r {params.Nirvana_ref} \
            --sd {params.Nirvana_supplementray} \
            -c {params.Nirvana_cache} {params.extra}
        """

rule Annovar:
    input:
        "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("05_Annotation/ANNOVAR/{type}")
    threads: config["annotation_threads"]
    resources:
        mem_mb=int(config["annotation_mem"])* 1024,
        cores=config["annotation_threads"],
        mem_gb=int(config["annotation_mem"]),
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
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
