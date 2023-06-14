
rule consequence:
    input: "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output: "04_calling/QC/{type}/bcftools_csq.vcf"
    threads: 1
    resources:
        mem_mb=2048,
        cores=1,
        mem_gb=2,
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

rule Nirvana:
    input:
        "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        "05_Annotation/Nirvana/{type}/Annotation.json.gz"
    threads: 1
    resources:
        mem_mb=8192,
        cores=1,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        Nirvana_supplementray = f"{nirvana_path}/DB/SupplementaryAnnotation/GRCh38/",
        Nirvana_ref = f"{nirvana_path}/DB/References/Homo_sapiens.GRCh38.Nirvana.dat",
        Nirvana_cache = f"{nirvana_path}/DB/Cache/GRCh38/Both",
        Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll",
        file_name = lambda wildcards: "05_Annotation/Nirvana/{wildcards.type}/Annotation"

    shell:
        """
        dotnet {params.Nirvana_cmd} \
            -i {input} \
            -o {params.file_name} \
            -r {params.Nirvana_ref} \
            --sd {params.Nirvana_supplementray} \
            -c {params.Nirvana_cache}
        """

rule Annovar:
    input:
        "04_calling/{type}/variants_genotyped.gvcf.gz"
    
    conda: "../env/wes_gatk.yml"

    output:
        directory("05_Annotation/ANNOVAR/{type}")
    threads: 4
    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        annovar_dir = annovar_dir,
        protocol = config["annovar_protocol"],
        operation = config["annovar_operation"],
        output = lambda wildcards: "05_Annotation/ANNOVAR/{wildcards.type}/annotations" 

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
            --thread {threads}
        """
