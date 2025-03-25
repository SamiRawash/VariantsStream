####installer

### download resource 
mkdir resource
cd resource
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz >  resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi > resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz > resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi > resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz.tbi
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta > resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta
curl https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai > resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai
curl https://raw.githubusercontent.com/hall-lab/speedseq/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed > exclude_hg38.bed
mkdir ./snpEff
mv snpEff
wget https://snpeff.blob.core.windows.net/databases/v5_0/snpEff_v5_0_GRCh38.99.zip
unzip snpEff_v5_0_GRCh38.99.zip
mv data
mkdir genomes
mv genomes
wget https://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz
mv Homo_sapiens.GRCh38.dna.toplevel.fa.gz GRCh38.99.fa.gz
cd ../../..
mkdir humandb
mv humandb

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    knownGene \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    refGene \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    cytoBand \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    avsnp150 \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    dbnsfp42c \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    gnomad312_genome \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    icgc28 \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    clinvar_20221231 \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    intervar_20180118 \
    ${PWD}/annovar/

singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    dbscsnv11 \
    ${PWD}/annovar/


singularity exec --bind "${PWD}" ../AnnotationFiltering/Annovar.sif \
    perl /usr/local/bin/annovar/annotate_variation.pl \
    -buildver hg38 \
    -downdb \
    -webfrom \
    annovar \
    dbnsfp31a_interpro \
    ${PWD}/annovar/
