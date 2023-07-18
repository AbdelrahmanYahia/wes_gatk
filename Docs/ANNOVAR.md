# ANNOVAR

- **Novagen main databases**
    - [ ]  RefSeq, Gencode → Gene
    - [ ]  All region filtration including: cytoband, snoRNA, miRNA, conservation regions, transcription factor binding site and repeats.
    - [ ]  1000 Human Genome, Exome Aggregation Consortium (ExAC) and exome sequencing project (ESP), were used to find alternative allele frequencies in populations that are reported.
    - [ ]  dbSNP
    - [ ]  COSMIC
    - [ ]  OMIM
    - [ ]  GWAS Catalog
    - [ ]  HGMD
    - [ ]  Gene Ontology
    - [ ]  KEGG
    - [ ]  Biocarta
    - [ ]  PID
    - [ ]  SIFT, PolyPhen, MutationAssessor, LRT and CADD scores were used to predict the deleteriousness of mutations. GERP++ scores were used to access the conservation of mutations.

- **Eqv. ANNOVAR database names**
    - [ ]  refGeneWithVer, knownGene, ensGene,
    - [ ]  phastConsElements46way, cytoBand, tfbsConsSites, wgRna, targetScanS, genomicSuperDups, dgvMerged, gwasCatalog, wgEncodeCaltechRnaSeqRawSignalRep1Gm12878CellLongpolyaBb12x75, wgEncodeRegTfbsClustered, snp153, wgEncodeUwDnaseSeqHotspotsRep2Gm12878
    - [ ]  1000g2015aug, exac10, esp6500siv2_all
    - [ ]  avsnp150
    - [ ]  COSMIC needs a custom script ([https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#cosmic-annotations](https://annovar.openbioinformatics.org/en/latest/user-guide/filter/#cosmic-annotations) )
    - [ ]  the OMIM data are the property of Johns Hopkins University and will not be available for download from UCSC (alterante link: [https://omim.org/downloads](https://omim.org/downloads) )
    - [ ]  GWAS
    - [ ]  No HGMD supplied by annovar
    - [ ]  No GO supplied by annovar
    - [ ]  No Kegg supplied by annovar
    - [ ]  No biocarta supplied by annovar
    - [ ]  PID ([http://www.openbioinformatics.org/annovar/download/GDI_full_10282015.txt.gz](http://www.openbioinformatics.org/annovar/download/GDI_full_10282015.txt.gz) )
    - [ ]  dbnsfp42c
- **Novagen annovar params:**
    
    This field tells whether the variant hits exons or hits intergenic regions, or hits introns, or hits non-coding RNA genes. The value of this field takes the following precedence: exonic = splicing > ncRNA > UTR5/UTR3 > intronic > upstream/downstream > intergenic. Notes: 1. When a variant hits different genes or transcripts, the variant may fit multiple functional categories, and then the precedence mentioned above is used to decide what function to print out; 2. The 'exonic' here refers only to coding exonic portion, but not UTR portion, as there are two keywords (UTR5, UTR3) that are specifically reserved for UTR annotations; 3. If a variant is located in both 5'UTR and 3'UTR region (possibly for two different genes), then the 'UTR5,UTR3' will be printed as the output; 4. 'splicing' in ANNOVAR is defined as variant that is within 2bp away from an exon/intron boundary by default; 5. 'splicing' in ANNOVAR only refers to the 2bp in the intron that is close to an exon; 6. The term 'upstream' and 'downstream' is defined as 1kb away from transcription start site or transcription end site, respectively, taking in account of the strand of the mRNA. If a variant is located in both downstream and upstream region (possibly for 2 different genes), then the 'upstream, downstream' will be printed as the output.
    

---

<aside>
➡️ It is important to explain the difference between region-based annotation and filter-based annotation here. Filter-based annotation looks exact matches between a query variant and a record in a database; two items are identical only if they have identical chromosome, start position, end position, ref allele and alaternative allele. Region-based annotation looks for over lap of a query variant with a region (this region could be a single position) in a database, and it does not care about exact match of positions, and it does not care about nucleotide identity at all. So in a sense, region-based annotation is somewhat similar to tabix, except that it does have involve index so it is much slower, yet it allows more user configuration to fine-tune results.

</aside>

### Params modifications:

- To use human genome variation society's format:
    
    ```bash
            mkdir -p outdir
            perl annovar_dir/table_annovar.pl input.vcf.gz annovar_dir/humandb/ \
                -buildver hg38 \
                -out outdir/outprefix -remove \
                -protocol refGene,avsnp150,clinvar_20221231 \
                -operation g,f,f \
    						-arg '-hgvs',,  \
                -nastring . \
                -vcfinput \
                --thread 12
    ```
    
    ```bash
    ~/annovar/table_annovar.pl variants_genotyped.filttered.gvcf.gz  ~/annovar/humandb/ -buildver hg38 -out test1 -protocol refGene,avsnp150,clinvar_20221231,cosmic70,dbnsfp31a_interpro,EAS.sites.2015_08,EUR.sites.2015_08,gme,gnomad211_exome,SAS.sites.2015_08,wgRna,refGeneVersion,gwasCatalog,cytoBand,ensGene,ensGeneMrna -operation g,f,f,f,f,f,f,f,f,f,r,g,r,r,g,r -arg '-hgvs -separate -exonsort',,,,,,,,,,,,,, -nastring . -vcfinput -remove -thread 12 -polish -xref ~/annovar/example/gene_fullxref.txt
    ```
    
    Here `-arg` takes a comma delimited list similar in size to both `-protocol` and `-operation` and each variable contains all the args for this operation like for `g` operation you can write  `'-hgvs -seperate'` . T*he `-transcript_function`* argument can be used to print transcript name instead of gene name.
    
- table_annovar automatically  converts the vcf file to anvin format using [`convert2annovar.pl](http://convert2annovar.pl/) -includeinfo -allsample -withfreq -format vcf4`
- in gene annotation to have the HGVS formatted strings for not only exonic variant, but also intronic variant that could be say 10bp away from splice site, by default, ANNOVAR only treats variants within 2bp of exon/intron boundary as splice variants, unless a --slicing_threshold parameter is set
- If you need the zygosity, quality and read coverage information in the output line as well, add the `-withzyg` argument with `convertarg`
- Gene definition: *1. If a gene is annotated as both coding and non-coding (multiple transcripts, some coding, some non-coding), the gene will be regarded as coding (the non-coding transcript definition will be ignored). 2. If a gene or a transcript has one or several non-coding definitions but without coding definition, it will be regarded as ncRNA in annotation output. 3. If a transcript maps to multiple locations as "coding transcripts", but some with complete ORF, some without complete ORF (that is, with premature stop codon), then the ones without complete ORF will be ignored. 4. If a transcript maps to multiple locations, all as "coding transcripts", but none has a complete ORF, then this transcript will not be used in exonic_variant_function annotation and the corresponding annotation will be marked as "UNKNOWN". 5. NEW in July 2014: If a transcript maps to multiple genomic locations, all mapping wil be used in the annotation process. Previously, only the "most likely" mapping will be used in annotation.*
- only exonic variants are annotated with gene function (nonsynonymous SNV, synonymous SNV, frameshift insertion, frameshift deletion, nonframeshift insertion, nonframeshift deletion, frameshift block substitution, nonframshift block substitution)
- When specifying amino acid changes, the specification always relates to a position for a transcript (not a "gene"). However, due to alternative splicing, if there are two or more transcripts that are all annotated for a gene, then the position of the amino acid change will differ, and it is important to always list the transcripts, in addition to gene names.
