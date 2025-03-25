Germline-NF PIPELINE
        Fabio Landuzzi, Agata Fant, Debora Charrance, Federica Furia e Stefano Marangoni, Progetto 5000genomi@VdA, CMP3@VdA (2023)


                                                                                       (((*                                                                             
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
                                                                                                                                                                


        5555555   000000    000000    000000                                                        ii      @@@@@@@      VV   VV         dd   AAAAAA 
        55       000   00  000   00  000   00                                                              @@@    @@@    VV   VV         dd  AA    AA
        55       0000  00  0000  00  0000  00  gggggg    eeeeee   nnnnnnn    oooooo   mmmmmm mmmm   ii   @@@        @@   VV   VV    ddddddd  AA    AA
        5555555  00 00 00  00 00 00  00 00 00 gg    gg  ee    ee  nn    nn  oo    oo  mm   mm   mm  ii  @@   @@@@@  @@   VV   VV   dd    dd  AAAAAAAA
              55 00  0000  00  0000  00  0000 gg    gg  eeeeeeee  nn    nn  oo    oo  mm   mm   mm  ii  @@  @@  @@  @@    VV VV    dd    dd  AA    AA
        55    55 00   000  00   000  00   000 gg    gg  ee        nn    nn  oo    oo  mm   mm   mm  ii  @@  @@  @@  @@     VVV     dd    dd  AA    AA
         555555   000000    000000    000000   ggggggg   eeeeeee  nn    nn   oooooo   mm   mm   mm  ii  @@   @@@@@@@@       V       ddddddd  AA    AA
                                                    gg                                                   @@                                          
                                               gg   gg                                                    @@@    @@@                                 
                                                gggggg                                                      @@@@@@                                  


        The typical command for running the pipeline is as follows:

        nextflow run Pipeline_Germline_DSL2.nf -entry {workflow} --samples samples_ids --outdir /path/to/results/ 

        PARAMETERS:
        	 -entry {workflow}: select the workflow of choice to be run, options are: 
                        PB_germ: Only Parabricks Germline, from fastq to bam and vcf 
                        Pipeline_Germline_SNV: Parabricks with SNV, from fastq to annotated+filtered SNV
                        Pipeline_Germline_SNV_nofilt: Parabricks with SNV, from fastq to annotated SNV
                        Pipeline_Germline_SNV_CNVkit: Pipeline with SNV and CNVkit, from fastq to annotated+filtered SNV and CNV from CNVkit  
                        Pipeline_Germline_SNV_CNVkit_SV: Pipeline with SNV, SV and CNVkit, from fastq to annotated+filtered SNV and CNV from CNVkit
                        CNVkit: CNVkit only, from bam
                        SNV_annot_filt: SNV annotation and filtration, from vcf to annotated+filtered SNV
                        SNV_annot: SNV annotation, from vcf to annotated SNV
                        SNV_filt: SNV filtration, from vcf to filtered SNV
                        SV_consensus: SV consensus pipeline, from bam and vcf to SV annotated
        	 
		 --samples: list of ids of the sample, comma separated, mandatory parameter

		 --outdir: path and directory name of the output directory, facultative parameter

HOW TO RUN THE PIPELINE AND HOW THE PIPELINE BEHAVE
	
	-Before starting the pipeline the user has to set the directory where the FASTQ file are present, in case of 5000genomi samples the location for all the projects are the same ("/archive/genomics_data/5000genomi/FASTQ/"),
	followed by "*/*_R{1,2}.fastq.gz". The first "*" take all the directory containing the fastq of the sample, the second "*" generalize the sample_id of the fastq.
	
	-To run the pipeline, the user can specify only the samples using the parameters "--samples", the lists has to be a comma separeted list of ids (ex. "ONC00001BLD01R,ONC0000FBLD01R").

	-the results of the pipeline are structered like this:
			+----> "sample" +----> "sample_id.bam/bam.bai" : bam from Parabricks
			|		|	
			|		+----> "sample_id.cnr/cns" : CNVkit results
			|		|
			|		+----> "sample_id.vcf" : VCF for Parabricks
			|		|
			|		+----> "sample_id_DuplMetrics.txt/_Recal.txt" : txt files containing the info for the Duplication metrics and recalibration
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
			 
