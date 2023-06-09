

rule consequence:
    input: "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output: "04_calling/QC/bcftools_csq.vcf"
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
        "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output:
        "05_Annotation/Nirvana/Annotation.json.gz"
    threads: 1
    resources:
        mem_mb=8192,
        cores=1,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        Nirvana_supplementray = f"{nirvana_path}/Data/SupplementaryAnnotation/GRCh38/",
        Nirvana_ref = f"{nirvana_path}/Data/References/Homo_sapiens.GRCh38.Nirvana.dat",
        Nirvana_cache = f"{nirvana_path}/Data/Cache/GRCh38/Both",
        Nirvana_cmd = f"{nirvana_path}/bin/Release/net*/Nirvana.dll",
        file_name = "05_Annotation/Nirvana/Annotation"

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
        "04_calling/variants_genotyped.gvcf.gz"
    
    conda: "GUAP"

    output:
        directory("05_Annotation/ANNOVAR")
    threads: 4
    resources:
        mem_mb=8192,
        cores=4,
        mem_gb=8,
        nodes = 1,
        time = lambda wildcards, attempt: 60 * 2 * attempt
    params:
        annovar_dir = "~/annovar",
        protocol = "refGene,avsnp150,clinvar_20221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,EUR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08",
        operation = "g,f,f,f,f,f,f,f,f,f",
        output = "05_Annotation/ANNOVAR/annotations"

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
