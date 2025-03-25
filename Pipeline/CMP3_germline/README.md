# Germline-NF PIPELINE
*Authors: Stefano Marangoni, Agata Fant, Debora Charrance, Federica Furia, Fabio Landuzzi e Alessandro Coppe Progetto 5000genomi@VdA, CMP3@VdA (2023)*

![logo](logo.png)
<!--                                                                                        (((*                                                                             
                                                                                    *((((((((.                                                                          
                                                                                  (((((((((((((*                                                                        
                                                                               *((((((|  *((((((((.                                                                     
                                                                             (((((((,       ((((((((*                                                                   
                                                                          *((((((|            *(((((((                                                                  
                                                          |(*            |(((((,                 (((((((,          .((.                                                 
                                                        .((((((       ,(((((|                      ,(((((((       ,(((((|                                               
                                                      *(((((((((,   |((((((*           (((*           (((((((,  |(((((((((                                              
                                                     |((((  ,((((((((((((.          *((((((((.          ,((((((((((|  *(((((*                                           
                                                 ,(((((|    .((((((((*           (((((    ((((*           (((((((,    ,(((((|                                          
                                                *((((.         *(((((*        *((((         (((((.        ((((((.        |(((((,                                       
                                             .((((((            ,(((((|    .(((((     ((((    (((((*    .((((|            *(((((*                                      
                             .,,,,,,,,,,,,,  ,(((((|              .((((((.*((((    ((((.*((((   ((((((((((((,              ,(((((|  .,,,,,,,,,,,,,                     
                       ,,,,,,,,,,,,,,,,,,,,,,,,**.                   *(((((((   ((((,        |((((((((((((((.                  **,,,,,,,,,,,,,,,,,,,,,,,,.              
                   .,,,,,,,,, ,(|     (((*  .,,,,,,,,                *(((((  (((*                (((((((((((.              .,,,,,,,, *(((.    |((| .,,,,,,,,,           
                 ,,,,,*|.     ,(|     (((*    .((*,,,,,.          .(((((((((|                      ,((((((((((|          ,,,,,*|.    *(((.    |((|     ((*,,,,,.        
              .,,,,,,*((.     ,(|     (((*    .(((, ,,,,,,      *(((((((((,                            *(((((((((.    .,,,,,*(((,    *(((.    ,***,    (((|,,,,,,,      
            ,,,,,,,  .((.   ,,,,,,    (((*    .(((,   .,,,,,..((((((((*                                    ((((((((|,,,,,,  ,(((,    *(((.   ,,,,,,    (((*  .,,,,,,.   
          .,,,,,.   .,,,,.   .,,,,    (((*     ,,,,.   .,,,,,*((((|                                          ,((((|,,,,,.  .,,,,,.   *(((.    .,,,    ,,,,,,    ,,,,,,  
        ,,,,,,     .,,,,,.   .,,,     ,**,,  .,,,,,.      ,,,,,,                                                ,,,,,,      .,,,,.  ,,,,,,    .,,,,    ,,,,,      .,,,,,
        ,,,,,,              ,,,,,,   ,,,,,,               ,,,,,,,.                                             ,,,,,,,      .,,,,.  .,,,,,   ,,,,,,    ,,,.       .,,,,,
         .,,,,,.    .,,,,.   .*|*,     ..    .,,,,,.     ,,,,,,,,,,,                                        .,,,,,,,,,,.   .,,,,,.            |((|    ,,,,,,     ,,,,,, 
            ,,,,,.  .,,,,.    ,(|     ,,,,,    ,,,,.  .,,,,,.  .,,,,,                                     ,,,,,,,   ,,,,,,  ,|||.   .,,,,,    |((|     |||*   .,,,,,.   
             ,,,,,,,..((.     ,(|    ,,,,,,     |(, ,,,,,,,       ,,,,,.                                .,,,,,.      ,,,,,,,*(((,    .,,,,    |((|     (((* ,,,,,,      
               .,,,,,*|(.     ,(|     (((*      |(*,,,,,.           .,,,,,,.                         ,,,,,,,           .,,,,,*((,    *(((.    |((|     ((|*,,,,,,       
                  ,,,,,,,,.   ,(|     (((*   .,,,,,,,,                 ,,,,,,,,.                 ,,,,,,,,.                ,,,,,,,,.  *(((.    |((|    ,,,,,,,,.         
                      ,,,,,,,,,**,,,,,***,,,,,,,,,.                      .,,,,,,,,,,,,,,,,,,,,,,,,,,,,                        ,,,,,,,,,**,,,,,,**,,,,,,,,,.             
                           ,,,,,,,,,,,,,,,,,,.                                .,,,,,,,,,,,,,,,,,,,                                 ,,,,,,,,,,,,,,,,,,.                  
                                                                                                                                                                


        5555555   00000    00000    00000                                                      ii      @@@@@@@     VV   VV        dd   AAAAA 
        55       000  00  000  00  000  00  gggggg    eeeeee  nnnnnnn    ooooo   mmmmm mmmmm        @@@       @@@  VV   VV        dd  AA   AA
        555555   00 0 00  00 0 00  00 0 00 gg    gg  ee    ee nn    nn  oo   oo  mm   mm   mm  ii  @@   @@@@@  @@  VV   VV    dddddd  AAAAAAA
             55  00  000  00  000  00  000 gg    gg  eeeeeee  nn    nn  oo   oo  mm   mm   mm  ii  @@  @@  @@  @@   VV VV   dd    dd  AA   AA
        55   55  00   00  00   00  00   00 gg    gg  ee       nn    nn  oo   oo  mm   mm   mm  ii  @@  @@  @@  @@    VVV    dd    dd  AA   AA
         55555    00000    00000    00000    gggggg   eeeeee  nn    nn   ooooo   mm   mm   mm  ii  @@   @@@@@@@@      V      ddddddd  AA   AA
                                                 gg                                                  @@@                                       
                                            gg   gg                                                    @@@@@@                               
                                             ggggg                                                                                        
 -->

ParGermPipe is a Pipeline which using Parabricks and Nextflow offers a highly Paralelized and efficient Pipeline for short-read technology for HPC ambient. Starting from the FASTQ files returns various type of annotated Variants (SNV, CNV and SV). The general workflow is mainly composed by 5 sub-workflows: 
- Clara Parabricks Germline: starting from the FASTQ files, perform the Allignment using BWA-MEM GPU port, and SNV calling using HaplotypeCaller GPU port, returning the BAM file and SNV VCF file.
- SNV Annotation: starting from SNV VCF file, perform the annotation of SNV variants using SNPEff, Annovar and Cosmic Databases, returning the Annotated SNV VCF file.
- SNV Filtration: starting from the Annotated SNV VCF file, perform a filtering based on the PASS flag, the Allele Frequency, the Impact factor added by SNPEff, the CancerGeneCensus filtering, a Intervar filtering and a Incidental Findings, returning the Annotated and Filtered SNV VCF file. 
- CNVkit: starting from the BAM file, perform CNV calling using CNVkit, returning the CNV VCF file and plot images. 
- SV: starting from the BAM file and SNV VCF file, perform SV calling using Manta, CNVnator, Lumpy and Breakdancer, the consensus calling of the before-mentioned caller and annotation with Duphold, AnnotSV and SVAFotate and filtering, returning the SV annotated and filtered VCF file. 


## Installation Guide

First of all, check that the dependecies are met: 
- Singularity
- make
- java, at least version 11

Than clone the repository.

      git clone https://gitlab.iit.it/Stefano.Marangoni/cmp3-germline-pipeline.git

After cloning the repository, run the Makefile installer.

      #if you want to install Nextflow locally, install the container and install the resources needed
      make 

      #if you want install Nexflow locally and the containers
      make container 

      #if you want to install the resources needed
      make download-database 

PS: Make sure to have the sudo privelege if you are installing the containers or ask to someone with sudo priveleg to install it for you.

## How to run ParGermPipe and how the pipeline behave

The typical command for running the pipeline is as follows:

      nextflow run Pipeline_Germline_DSL2.nf -entry {workflow} --samples samples_ids --outdir /path/to/results/ 

PARAMETERS:

- -entry {workflow}: select the workflow of choice to be run, facultative parameter (default= None, if not specified it will run the pipeline with SNV, SV and CNVkit, from fastq to annotated+filtered SNV, SV and CNV from CNVkit) 
    - PB_germ: Only Parabricks Germline, from fastq to bam and vcf 

    - Pipeline_Germline_SNV: Parabricks with SNV, from fastq to annotated+filtered SNV
    
    - Pipeline_Germline_SNV_nofilt: Parabricks with SNV, from fastq to annotated SNV
    
    - Pipeline_Germline_SNV_CNVkit: Pipeline with SNV and CNVkit, from fastq to annotated+filtered SNV and CNV from CNVkit  
    
    - CNVkit_workflow: CNVkit only, from bam
    
    - Pipeline_Germline_SNV_nofilt_SV: Pipeline with SNV and SV, from fastq to annotated SNV and annotated and filtered SV
    
    - Pipeline_Germline_SNV_nofilt_CNVkit_SV: Pipeline with SNV and SV, from fastq to annotated SNV, annotated and filtered SV Pipeline with SNV and SV, from fastq to annotated SNV and SV 
		
    - Pipeline_Germline_SV:  Pipeline with SNV and SV, from fastq to SNV and SV
	  
    - Pipeline_Germline_CNVkit_SV:  Pipeline with SNV, SV and CNVkit, from fastq to SNV and CNV from CNVkit
		
    - Coverage_workflow: Coverage process only, from bam

    - AddPassHardFilter: AddPass process only, from vcf without PASS flag to vcf with PASS flag

    - SNV_AddPass_annot_filt:

    - SNV_annot_filt: SNV annotation and    filtration, from vcf to annotated and filtered SNV
    
    - SNV_annot: SNV annotation, from vcf to annotated SNV 
    
    - SNV_filt: SNV filtration, from vcf to filtered SNV
    
    - SV_consensus: SV consensus pipeline, from bam and vcf to SV annotated and filtrated 
		
    - SV_consensus_no_MELT:  SV consensus pipeline, from bam and vcf to SV annotated and filtrated without MELT

    - SV_annot_filt: SV annotation and Filtration, from SV vcf to SV annotated and filtrated

    - SV_annot_filt_no_Melt: SV annotation and Filtration, from SV vcf to SV annotated and filtrated without MELT

    - Melt_workflow: MELT process only, from bam


- --samples: list of ids of the sample, comma separated, mandatory parameter

- --outdir: path and directory name of the output directory, facultative parameter

Before starting any pipeline the user has to set the directory where the FASTQ/BAM/VCF files are present ("/path/to/location_of/FASTQ/") followed by "\*/\*_R{1,2}.fastq.gz". The first "\*" take all the directory containing the fastq of the sample, the second "\*" generalize the sample_id of the fastq.
  
- Decide the pipeline to run selecting it with the parameter "-entry", the possible pipeline runnable are listed above.
  
- To run the pipeline, the user can specify only the samples using the parameters "--samples", the lists has to be a comma separeted list of ids (ex. "ONC00001BLD01R,ONC0000FBLD01R").
  
- Based on the required input of the pipeline, the FASTQ/BAM/VCF are taken from the path given automatically

- The results of the pipeline are structered like this:


      +----> "sample" +----> "sample_id.bam/bam.bai" : bam from Parabricks    
      |		            |	         
      |		            +----> "sample_id.cnr/cns" : CNVkit results         
      |		            |           
      |		            +----> "sample_id.vcf" : VCF for Parabricks     
      |		            |                                                                           
      |		            +----> "sample_id_DuplMetrics.txt/_Recal.txt" : txt files containing the info for the Duplication metrics and recalibration     
      |       
      +----> "00_Coverage"        
      |        
      +----> "01_SnpEff" +----> "sample_id_SnpEff.csv/.vcf"	            	
      |        
      +----> "02_cosmic_annotation" +----> "sample_id_COSMIC.vcf"      
      |       
      +----> "03_filter_af" +----> "sample_id_AF_filtered.vcf"      
      |      
      +----> "04_filter_impact" +----> "sample_id_impact_filtered.vcf"      
      |      
      +----> "05_filter_pass" +----> "sample_id_PASS_filterd.vcf"     
      |      
      +----> "06_cancer_gene_census" +----> "sample_id_CGC.vcf"     
      |    
      +----> "07_InterVar": not yet implemented, soon implemented    
      |     
      +----> "08_TablaofGermline": not yet implemented, soon implemented     
      |     
      +----> "09_IncidentalFindings_Tables" +----> "sample_id.tsv"     
      |     
      +----> "Manta" +----> "sample_id_mantaGerm_INVconverted.vcf"    
      |     
      +----> "Breakdancer" +----> "sample_id_S4.vcf"     
      |     
      +----> "CNVnator" +----> "sample_id_S10.vcf"     
      |     
      +----> "Lumpy" +----> "sample_id-smoove.genotyped.vcf.gz"      
      |     
      +----> "SV" +----> "sample_id-smoove.genotyped_Duphold.vcf": Duphold output vcf file     
      |     
      +----> "sample_id_consensusFiltered_final.vcf": Duphold filtered output vcf file     
      |     
      +----> "sample_id.tsv" : AnnotSV annoted variants   
      |    
      +----> "sample_id.unannotated.tsv" : AnnotSV unannoted variants   

## PARAMETERS LIST
### PARAMETER FOR COMMON MANAGEMENT
This parameters are common to be used and can be found in the *nextflow.config* file. These parameters can be modified from 'nextflow.config' or as a parameters in the command line. In order to pass the parameter thourgh command line the path has to be included in double quotes ("").

- --outdir = path where the links results files are placed, by default the path is the the path where the Pipeline is launched

- --fastq_path = path of the fastq files for the samples, add "*\*/\*_R{1,2}_001.fastq.gz*" in order to make the pipeline gather from all the directory all the fastq pairs. 

- --bam_sample = path of the bam files for the samples, add "*\*/\*.bam{,.bai}*" in order to make the pipeline gather from all the directory all the bam/bai pairs.    

- --vcf_path = path of the SNV vcf files for the samples, add "*\*/\*_PassHF.vcf.gz{,.tbi}*" in order to make the pipeline gather from all the directory all the bam/bai pairs.

- --snv_annot_path = path of the annotated SNV vcf files for the samples, add "*\*_PassHF_MultiallelicSplit_SnpEff_annovar_COSMIC.vcf*" in order to make the pipeline gather from all the directory all the bam/bai pairs.    

- --sv_path = path of the SV vcf files for the samples, add "*\*_consensusFiltered_final.vcf*" in order to make the pipeline gather from all the directory all the bam/bai pairs.    


In the *nextflow.config* file there are also other parameters that are less common to be modified: 

- --modules_location = path where all the file with the processes are stored, by default "${baseDir}/bin/germ_process"

- --container_dir = path where all the containers for the processes are stored, by default "${baseDir}/bin/def"

- --resource_dir = path where all the recourses for the processes are stored such as refereces, annotation databases and others, by default "${baseDir}/bin/resource"

- --script_dir= path where all the scripts used in  the processes are stored, by default "${baseDir}/bin/script"

### PARAMETER FOR SPECIFIC MANAGEMENT
These parameters are for specify the container paths, usually contained in the "${baseDir}/bin/def" path (specified in the parameter *--container_dir* ).

- --PARABRICKS_4 = path to container for Parabricks 4, by default *\${params.container_dir}/parabricks_4.1/clara-parabricks_4.sif* 

- --CONTAINER_GATK = path to container for GATK, by default *\${params.container_dir}/parabricks_4.1/clara-parabricks_4.sif*  

- --CONTAINER_BCFTOOLS = path to container for bcftools, by default *\${params.container_dir}/parabricks_4.1/clara-parabricks_4.sif*   

- --CONTAINER_ANNOVAR = path to container for Annovar, by default *\${params.container_dir}/AnnotationFiltering/Annovar.sif*     

- --CONTAINER_SNPEFF = path to container for SNPEff, by default *\${params.container_dir}/AnnotationFiltering/SnpEff5_1.sif*

- --CONTAINER_FILTERING_AND_TABLES = path to container for tabularization and filtering, by default *\${params.container_dir}/VCF_AnnotationFiltering/VCF_AnnotationFiltering_v9.sif*     

- --CONTAINER_CNVKIT = path to container for CNVKit, by default *\${params.container_dir}/cnvkit/cnvkit.sif*      

- --CONTAINER_CNVKIT_R = path to container for CNVKit plot, by default *\${params.container_dir}/cnvkit/r_with_ps.sif*       

- --CONTAINER_MANTA = path to container for Manta, by default *\${params.container_dir}/manta/strelkaandmanta.sif*      

- --CONTAINER_LUMPY = path to container for Lumpy, by default *\${params.container_dir}/smoove/smoove_latest.sif*      

- --CONTAINER_SVE = path to container for SVE, by default *\${params.container_dir}/sve/sve_latest.sif*       

- --CONTAINER_SURVIVOR = path to container for Survivor, by default *\${params.container_dir}/survivor/survivor.sif*         

- --CONTAINER_SV_GENERAL = path to container for SV Pipeline various minor processes, by default *\${params.container_dir}/sv_cons_general_utilities/sv_cons_general_utilities.sif*    

- --CONTAINER_ANNOTSV = path to container for AnnotSV, by default *\${params.container_dir}/annotSV/annotSV_3_2_3.sif*   

- --COVERAGE_CONTAINER = path to container for Samplot used to calcolate the coverage of the bam file, by default *\${params.container_dir}/SamBcfBed_Cov/SamBcfBed_Tools.sif* 

- --SVAFotate_container = path to container for, by default *\${params.container_dir}/SVAFotate/SVAFotate.sif*    

- --Filtering_container = path to container for, by default *\${params.container_dir}Filtering/TSV_FilterAnnotation_MantaMelt_AnnotSV3.2.3_SVAfotate.sif* 



These other parameters are for resources needed by the processes:

- --ex_bed_smoove = path to bed file list of excluded sites used by smoove, by default *\${params.resource_dir}/exclude_hg38.bed*

- --ReferenceGenome = location of the Reference genome, by default *\${params.resource_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta*    

- --KnowSites = location of the Reference known sites of SNPs, by default *\${params.resource_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz* 

- --REF_HAPMAP_FILTERGERMLINE =  location of the resources of the hapmap, by default *\${params.resource_dir}/resources_broad_hg38_v0_hapmap_3.3.hg38.vcf.gz*

- --REF_MILLS_FILTERGERMLINE = location of the resources of the Mills, by default *\${params.resource_dir}/resources_broad_hg38_v0_Mills_and_1000G_gold_standard.indels.hg38.vcf.gz*

- --SNPEFF_DATABASES = location of the SNPEff databases, by default *\${params.resource_dir}/snpEff/data/*  

- --SNPEFF_REFERENCE_GENOME = filename of the SNPEff reference genome, by default *GRCh38.99*  

- --ANNOVAR_DATABASES =  location of the Annovar databases, by default *\${params.resource_dir}/annovar/humandb/*

- --ANNOVAR_REFERENCE_GENOME = filename of the Annovar reference genome, by default *hg38*

- --ANNOVAR_BED_FILE = filename of the Annovar bed file, by default *lncipedia_5_2_hg38.bed*

- --ANNOVAR_OMIM =  filename of the Annovar OMIM file, by default *genemap2_onlyGenes*

- --COSMIC_ANN_VCF = location of the resource for COSMIC annotation, by default *\${params.resource_dir}/cosmic/CosmicCodingMuts.vcf.gz*

- --CGC_GENE_FILE = location of the resource for COSMIC annotation, by default *"\${params.resource_dir}/cosmic/final_PanelApp_CGC_genes.txt"*

- --BED_SVAFotate = location of the resource for SVAFotate annotation, by default *${params.resource_dir}/SVAFotate_core_SV_popAFs.GRCh38.bed*

- --BED_for_filtering = location of the resource for filtering, by default *${params.resource_dir}/final_PanelApp_CGC_genes.bed*

- --CNVKIT_COVERAGE_REFERENCE = location of the reference for CNVkit, by default *${params.resource_dir}/reference.cnn*   

- --CNVKIT_centromere_position = location of the reference for CNVkit, by default *${params.resource_dir}/centromere_position.tsv*

- --CNVKIT_gene_to_filter = location of the list of gene, by defualt *${params.resource_dir}/cancer_gene_list.txt*

- --GENELIST_AND_CANONICALTR = location of the list of gene for the Table of Germline process, by defualt "${params.resource_dir}/PanelApp_CGC_genes.txt"



These parameters are used to assign the name of the folder of the results:  

- --COVERAGE_FOLDER = name of the directory containin the output file from Coverage process, by default *"00_Coverage"*

- --SNPEFF_FOLDER = name of the directory containin the output file from SNPEff process, by default *"01_SnpEff"*

- --COSMIC_FOLDER = name of the directory containin the output file from Annovar process, by default *"02_cosmic_annotation"*

- --FILTER_AF_FOLDER = name of the directory containin the output file from filtering based on AF process, by default *"03_filter_af"*

- --FILTER_IMPACT_FOLDER = name of the directory containin the output file from filtering based on Impact process, by default *"04_filter_impact"*

- --FILTER_PASS_FOLDER = name of the directory containin the output file from filtering based on PASS flag process, by default *"05_filter_pass"*

- --CGC_FOLDER = name of the directory containin the output file from filtering based on Cancer Gene Census process, by default *"06_cancer_gene_census"*

- --INTERVAR_FOLDER = name of the directory containin the output file from filtering based on Intervar process, by default *"07_InterVar"*

- --TABLES_GERM_FOLDER = name of the directory containin the output file from Table of Germline process, by default *"08_Tables_Germline"*

- --TABLES_IncidentalFindings_FOLDER = name of the directory containin the output file from Table of Incidental Findings process, by default *"09_IncidentalFindings_Tables"* 

- --SVAFOTATE_FOLDER = name of the directory containin the output file from SVAFotate process, by default *"SVAFotate"*     

- --ANNOTSV_FOLDER = name of the directory containin the output file from AnnotSV process, by default *"AnnotSV"*      

- --TSVFILTERING_FOLDER = name of the directory containin the output file from filtering on TSV process, by default *"FilteredTSV"*  

- --dupl_metrics_folder= name of the directory containin the output file from MultiQC process, by default *"DuplMetrics"*  


These parameters are used to give the paths of scripts used in processes:

- --script = location of the script for Coverage process , by default *"\${params.script_dir}/BASH_PBSSamtoolsCoverage.bash"*

- --SNPEFF_EXEC =  location of the executable for SNPEff process , by default *"/usr/local/bin/snpEff/snpEff.jar"*   

- --ANNOVAR_EXEC = location of the executable for Annovar process , by default *"/usr/local/bin/annovar/table_annovar.pl"*   

- --SNPSIFT_PATH = location of the executable of SnpSift for COSMIC, AF filtering and Impact Filtering, by default *"/usr/local/bin/snpEff/SnpSift.jar"*

- --GET_PASS_VARIANTS_PATH = location of the executable for the PASS flag filtering process, by default *"\${params.container_dir}/AnnotationFiltering/scripts/filter_pass/get_pass_variants"*

- --CGC_PYTHON_SCRIPT_PATH = location of the script for Cancer Gene Census filtering process , by default *"\${params.container_dir}/AnnotationFiltering/scripts/cancer_gene_census_filter_python_script/cancer_gene_census_filter.py"*

- --INTERVAR_PY_SCRIPT_PATH = location of the executable of SnpSift for Intervar filtering process, by default *"\${params.script_dir}/InterVar_Clinvar_CADD_filtering.py"* 

- --CNVKIT_EXEC = location of the executable for CNVkit process, by default *"/usr/local/bin/cnvkit.py"* 

- --CNVKIT_FILTERING_SCRIPT = location of the executable for CNVkit filtering process, by default *"/opt/my_python/filter_log2_and_genes_from_cns.py"*   

- --CNVKIT_FILTERING_MODIFY_SCRIPT = location of the executable for modify CNVkit filtering process, by default *"/opt/my_python/modify_cnr.py"*

- --CNVKIT_PLOT_RSCRIPT = location of the executable for CNVkit plot process, by default *"\${params.script_dir}/CNVKIT_graph.R"*


These parameters are used to assign the computational resources for some processes:


- --SNPEFF_MEMORY =  quantity of memory given to the SNPEff process, by default *"7"*    

- --ANNOVAR_CPUS = quantity of processing thread given to the Annovar process, by default *"4"* 

- --ANNOVAR_MEMORY = quantity of memory given to the Annovar process, by default *"20"*

- --SNPSIFT_MEMORY = quantity of memory given to SNPSift for COSMIC, AF filtering and Impact Filtering process, by default *"4"*    

- --CNVKIT_CPUS = quantity of processing thread given to the CNVkit process, by default *"10"* 

- --CNVKIT_MEMORY =  quantity of memory given to the CNVkit process, by default *"30"*


These parameters are used to assign values for certain process:

- --AF =  type of field given for the AF filtering, by default *"AF"*

- --AF_VALUE =  value of Allele Frequency given as threshlod for the AF filtering, by default *"0.05"*

- --CGC_SEPARATOR = type of separator given for the Cancer Gene Census process, by default *","*    

- --CGC_COLUMN = number of the column given for the Cancer Gene Census process, by default *"1"*

- --CGC_ANNOTATOR = type of value annotator given for the Cancer Gene Census process, by default *"s"*

- --CGC_EXCLUDE =   type of value of exclude given for the Cancer Gene Census process, by default *""*

- --table_fields =  fields to be given for to the Table of Germline process, by default *"'AF,CADD_phred,InterVar_automated,CLNSIG,CLNREVSTAT,Func.refGene_latest,ExonicFunc.refGene_latest'"*

- --sequencing_lib =  value of the library used to be given to the Parabricks germline pipeline process, by default *"PCR-Free"* 

- --sequencing_platform =  value of the platform to be given to the Parabricks germline pipeline process, by default *"ILLUMINA"*