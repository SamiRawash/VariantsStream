#! /hpcshare/genomics/ASL_ONC/NextFlow_RunningDir/nextflow-22.04.4-all
nextflow.enable.dsl=2

def helpMessage() {
    log.info'''
        Germline-NF PIPELINE
        Stefano Marangoni, Agata Fant, Debora Charrance, Federica Furia, Fabio Landuzzi e Alessandro Coppe, Progetto 5000genomi@VdA, CMP3@VdA (2023)


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
                
                -entry {workflow}: select the workflow of choice to be run, facultative parameter (default= None, if not specified it will run the pipeline with SNV, SV and CNVkit, from fastq to annotated+filtered SNV, SV and CNV from CNVkit) 
                    options are: 
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

                        - SNV_AddPass_annot_filt: AddPass SNV, SNV annotation and filtration, from vcf to annotated and filtered SNV

                        - SNV_annot_filt: SNV annotation and filtration, from vcf to annotated and filtered SNV
                        
                        - SNV_annot: SNV annotation, from vcf to annotated SNV 
                        
                        - SNV_filt: SNV filtration, from vcf to filtered SNV
                        
                        - SV_consensus: SV consensus pipeline, from bam and vcf to SV annotated and filtrated 
                            
                        - SV_consensus_no_MELT:  SV consensus pipeline, from bam and vcf to SV annotated and filtrated without MELT

                        - SV_annot_filt: SV annotation and Filtration, from SV vcf to SV annotated and filtrated

                        - SV_annot_filt_no_Melt: SV annotation and Filtration, from SV vcf to SV annotated and filtrated without MELT

                        - Melt_workflow: MELT process only, from bam

                --samples: list of ids of the sample, comma separated, mandatory parameter

                --outdir: path and directory name of the output directory, default value launchDir, facultative parameter
'''.stripIndent()    
}

if (params.help) {
    helpMessage()
    exit 0
}

include {PipelineGermline; PipelineGermline_gvcf} from "${params.modules_location}/parabricks_germline.nf"
include {Coverage} from "${params.modules_location}/coverage.nf"
include {CopyFollowingLinks; MultiQC_DuplMetr} from "${params.modules_location}/duplmetr.nf"
include {CNVkit; CNVkit_filter; CNVkit_plot} from "${params.modules_location}/cnvkit.nf"
include {AddPass_VariantFiltration; AddPass_VariantFiltration_gz} from "${params.modules_location}/add_pass.nf"
include {MultiAllelicSplit} from "${params.modules_location}/multiallelic_split.nf"
include {SnpEff_Annotation; MultiQC_SnpEff} from "${params.modules_location}/snpeff.nf"
include {Annovar_Annotation} from "${params.modules_location}/annovar.nf"
include {CosmicAnnotation} from "${params.modules_location}/cosmic.nf"
include {IncidentalFindings_filter; Vcf2Table_IF} from "${params.modules_location}/incidental_fidings.nf"
include {AF_Filter} from "${params.modules_location}/af_filt.nf"
include {Impact_Filter} from "${params.modules_location}/impact_filt.nf"
include {PASS_Filter} from "${params.modules_location}/pass_filt.nf"
include {CGC_Filter} from "${params.modules_location}/cgc_filt.nf"
include {InterVar_Filter} from "${params.modules_location}/intervar.nf"
include {Vcf2Table} from "${params.modules_location}/table_germ.nf"
include {Manta} from "${params.modules_location}/manta.nf"
include {Lumpy} from "${params.modules_location}/lumpy.nf"
include {Cnvnator} from "${params.modules_location}/cnvnator.nf"
include {Breakdancer} from "${params.modules_location}/breakdancer.nf"
include {VCF_format_mod} from "${params.modules_location}/format_mod.nf"
include {Consensus_preprocess; Manta_ins_vcf} from "${params.modules_location}/del_dup_inv_tra_ins_preproc.nf"
include {Survivor_del_dup; Survivor_inv_tra; Fix_survivor_call} from "${params.modules_location}/survivor.nf"
include {Vcf_germ_ind; Duphold_no_tra; Postprocess_duphold} from "${params.modules_location}/duphold.nf"
include {AnnotSV} from "${params.modules_location}/annotsv.nf"
include {SVAFotate} from "${params.modules_location}/svafotate.nf"
include {TSV_filtering; TSV_filtering_no_MELT} from "${params.modules_location}/sv_filt.nf"
include {Melt; SVAFotate_Melt; AnnotSV_Melt; Melt_concatenate} from "${params.modules_location}/melt.nf"

lista=params.samples?.tokenize(',')

def sample_id_position(ids, auto_channel_id) {
        split_id=auto_channel_id.split('_')
        //println(split_id)
        selected_ind=0
        split_id.eachWithIndex{ id, ind -> if(ids.contains(id)){selected_ind=ind}}

        return selected_ind
}

fastq_ch=Channel.fromFilePairs(params.fastq_path, flat: true).filter{lista.contains(it[0].split('_')[sample_id_position(lista,it[0])])}
bam_ch=Channel.fromFilePairs(params.bam_path,flat: true, size:-1).filter{lista.contains(it[0].split('_')[sample_id_position(lista,it[0])])}.ifEmpty("No_Bam")
germ_snv_ch=Channel.fromFilePairs(params.vcf_path, flat: true, size:-1).filter{lista.contains(it[0].split('_')[sample_id_position(lista,it[0])])}.ifEmpty(["No Germ",["$baseDir/nextflow.config"]])
annot_ch=Channel.fromFilePairs(params.snv_annot_path, flat: true, size:-1).filter{lista.contains(it[0].split('_')[sample_id_position(lista,it[0])])}.ifEmpty(["No Annot",["$baseDir/nextflow.config"]])
germ_sv_ch=Channel.fromFilePairs(params.sv_path, flat: true, size:-1).filter{lista.contains(it[0].split('_')[sample_id_position(lista,it[0])])}.ifEmpty("No SV")

workflow pb_germline {
    take: fq_ch
    main:
        parabricks_results=PipelineGermline(fq_ch, file("${params.ReferenceGenome}"), file("${params.KnowSites}"),"${params.outdir}")
    
        bam_ch=parabricks_results[0]
        vcf_ch=parabricks_results[1]

        dir_DuplsMetr=file("${params.outdir}/${params.dupl_metrics_folder}")
        if ( ! dir_DuplsMetr.isDirectory() ) {
            result = dir_DuplsMetr.mkdir()
            println result ? "${dir_DuplsMetr}: OK" : "Cannot create directory: $dir_DuplsMetr"
        }

        CopyFollowingLinks(parabricks_results[3],file("${dir_DuplsMetr}"))
        MultiQC_DuplMetr(CopyFollowingLinks.out.collect(), "${dir_DuplsMetr}/${params.dupl_metrics_folder}")
    emit:
        bam_ch
        vcf_ch
}

workflow pb_germline_gvcf {
    take: fq_ch
    main:
        parabricks_results=PipelineGermline_gvcf(fq_ch, file("${params.ReferenceGenome}"), file("${params.KnowSites}"),"${params.outdir}")

        bam_ch=parabricks_results[0]
        vcf_ch=parabricks_results[1]

        dir_DuplsMetr=file("${params.outdir}/${params.dupl_metrics_folder}")
        if ( ! dir_DuplsMetr.isDirectory() ) {
            result = dir_DuplsMetr.mkdir()
            println result ? "${dir_DuplsMetr}: OK" : "Cannot create directory: $dir_DuplsMetr"
        }

        CopyFollowingLinks(parabricks_results[3],file("${dir_DuplsMetr}"))
        MultiQC_DuplMetr(CopyFollowingLinks.out.collect(), "${dir_DuplsMetr}/${params.dupl_metrics_folder}")
    emit:
        bam_ch
        vcf_ch
}

workflow coverage {
    take: bam_ch
    main:
        dir_Coverage=file("${params.outdir}/${params.COVERAGE_FOLDER}")
        if ( ! dir_Coverage.isDirectory() ) {
            result = dir_Coverage.mkdir()
            println result ? "${dir_Coverage}: OK" : "Cannot create directory: $dir_Coverage"
        }

        coverage_ch= Coverage(bam_ch)
    emit:
        coverage_ch
}

workflow add_pass{
    take: vcf_ch
    main:
	ind_ch=vcf_ch.branch{
                gz:     it.size()==3
                no_ind: it.size()==2
        }

    add_pass_gz_ch=AddPass_VariantFiltration_gz(ind_ch.gz, "2.0", "30.0", "3.0", "60.0", "40.0", "-12.5", "-8.0", "${params.outdir}")

    add_pass_nogz_ch=AddPass_VariantFiltration(ind_ch.no_ind, "2.0", "30.0", "3.0", "60.0", "40.0", "-12.5", "-8.0", "${params.outdir}")

    add_pass_ch=add_pass_nogz_ch.mix(add_pass_gz_ch).ifEmpty(add_pass_gz_ch.mix(add_pass_nogz_ch))

    emit:
        add_pass_ch
}

workflow germline_snv_annotation{
    take:
        add_pass_ch
    main:
        mas_ch=MultiAllelicSplit(add_pass_ch)

        dir_SnpEff=file("${params.outdir}/${params.SNPEFF_FOLDER}")
        if ( ! dir_SnpEff.isDirectory() ) {
            result = dir_SnpEff.mkdir()
            println result ? "${dir_SnpEff}: OK" : "Cannot create directory: $dir_SnpEff"
        }

        snpeff_ch=SnpEff_Annotation(mas_ch, "${params.SNPEFF_REFERENCE_GENOME}")
        vcf_snpeff_ch=snpeff_ch[0]
        MultiQC_SnpEff(snpeff_ch[1].collect())

        annovar_ch=Annovar_Annotation(vcf_snpeff_ch, "${params.ANNOVAR_REFERENCE_GENOME}", "${params.ANNOVAR_BED_FILE}","${params.ANNOVAR_OMIM}")


        if (params.COSMIC_FLAG) {
            dir_cosmic=file("${params.outdir}/${params.COSMIC_FOLDER}")
            if( ! dir_cosmic.isDirectory() ) {
                result = dir_cosmic.mkdir()
                println result ? "${dir_cosmic}: OK" : "Cannot create directory: $dir_cosmic"
            }

            cosmic_vcf_ch=CosmicAnnotation(annovar_ch[0],"${params.COSMIC_ANN_VCF}","${params.COSMIC_ANN_VCF}.tbi", "${params.SNPEFF_REFERENCE_GENOME}")
        } else{
            cosmic_vcf_ch=annovar_ch
        }
    emit:
        cosmic_vcf_ch
}

workflow cnvkit{
    take:
        bam_ch
    main:
        cnvkit_ch=CNVkit(bam_ch,"${params.CNVKIT_COVERAGE_REFERENCE}", "${params.CNVKIT_EXEC}","${params.outdir}")
        CNVkit_filter(cnvkit_ch, "${params.CNVKIT_gene_to_filter}", "${params.CNVKIT_FILTERING_SCRIPT}", "${params.CNVKIT_FILTERING_MODIFY_SCRIPT}","${params.outdir}")
        CNVkit_plot(CNVkit_filter.out,"${params.CNVKIT_PLOT_RSCRIPT}","${params.CNVKIT_centromere_position}","${params.outdir}")
        CNVkit_plot.out.count().map { print "CNVkit Plot Generated: $it" }
    emit:
	cnvkit_ch
}

workflow snv_filt{
    take: vcf_cosmic_ch
    main:
        dir_AF_Filter=file("${params.outdir}/${params.FILTER_AF_FOLDER}")
        if ( ! dir_AF_Filter.isDirectory() ) {
            result = dir_AF_Filter.mkdir()
            println result ? "${dir_AF_Filter}: OK" : "Cannot create directory: $dir_AF_Filter"
        }
        
        af_ch=AF_Filter(vcf_cosmic_ch, "${params.AF}","${params.AF_VALUE}")

        dir_Impact_Filter=file("${params.outdir}/${params.FILTER_IMPACT_FOLDER}")
        if ( ! dir_Impact_Filter.isDirectory() ) {
            result = dir_Impact_Filter.mkdir()
            println result ? "${dir_Impact_Filter}: OK" : "Cannot create directory: $dir_Impact_Filter"
        }
        
        impact_ch=Impact_Filter(af_ch)

        dir_PASS_Filter=file("${params.outdir}/${params.FILTER_PASS_FOLDER}")
        if ( ! dir_PASS_Filter.isDirectory() ) {
            result = dir_PASS_Filter.mkdir()
            println result ? "${dir_PASS_Filter}: OK" : "Cannot create directory: $dir_PASS_Filter"
        }

        pass_ch=PASS_Filter(impact_ch,"${params.GET_PASS_VARIANTS_PATH}")

        dir_CGC_Filter=file("${params.outdir}/${params.CGC_FOLDER}")
        if ( ! dir_CGC_Filter.isDirectory() ) {
            result = dir_CGC_Filter.mkdir()
            println result ? "${dir_CGC_Filter}: OK" : "Cannot create directory: $dir_CGC_Filter"
        }
        
        cgc_ch=CGC_Filter(pass_ch, "${params.CGC_PYTHON_SCRIPT_PATH}","${params.CGC_GENE_FILE}")
        
        dir_InterVar_Filter=file("${params.outdir}/${params.INTERVAR_FOLDER}")
        if ( ! dir_InterVar_Filter.isDirectory() ) {
            result = dir_InterVar_Filter.mkdir()
            println result ? "${dir_InterVar_Filter}: OK" : "Cannot create directory: $dir_InterVar_Filter"
        }

        intervar_ch=InterVar_Filter(cgc_ch,"${params.INTERVAR_PY_SCRIPT_PATH}")

        dir_Tables_Germ=file("${params.outdir}/${params.TABLES_GERM_FOLDER}")
        if ( ! dir_Tables_Germ.isDirectory() ) {
            result = dir_Tables_Germ.mkdir()
            println result ? "${dir_Tables_Germ}: OK" : "Cannot create directory: $dir_Tables_Germ"
        }

        Vcf2Table(intervar_ch, file("${params.GENELIST_AND_CANONICALTR}"))

        dir_Tables_IF=file("${params.outdir}/${params.TABLES_IncidentalFindings_FOLDER}")
        if ( ! dir_Tables_IF.isDirectory() ) {
            result = dir_Tables_IF.mkdir()
            println result ? "${dir_Tables_IF}: OK" : "Cannot create directory: $dir_Tables_IF"
        }

        if_ch=IncidentalFindings_filter(vcf_cosmic_ch)
        Vcf2Table_IF(if_ch,"${params.GENELIST_IF}")
}

workflow snv_filt_no_onc{
    take: vcf_cosmic_ch
    main:
        dir_AF_Filter=file("${params.outdir}/${params.FILTER_AF_FOLDER}")
        if ( ! dir_AF_Filter.isDirectory() ) {
            result = dir_AF_Filter.mkdir()
            println result ? "${dir_AF_Filter}: OK" : "Cannot create directory: $dir_AF_Filter"
        }
        
        af_ch=AF_Filter(vcf_cosmic_ch, "${params.SNPSIFT_PATH}", "${params.AF}","${params.AF_VALUE}", "${dir_AF_Filter}")

        dir_Impact_Filter=file("${params.outdir}/${params.FILTER_IMPACT_FOLDER}")
        if ( ! dir_Impact_Filter.isDirectory() ) {
            result = dir_Impact_Filter.mkdir()
            println result ? "${dir_Impact_Filter}: OK" : "Cannot create directory: $dir_Impact_Filter"
        }
        
        impact_ch=Impact_Filter(af_ch, "${params.SNPSIFT_PATH}")

        dir_PASS_Filter=file("${params.outdir}/${params.FILTER_PASS_FOLDER}")
        if ( ! dir_PASS_Filter.isDirectory() ) {
            result = dir_PASS_Filter.mkdir()
            println result ? "${dir_PASS_Filter}: OK" : "Cannot create directory: $dir_PASS_Filter"
        }

        pass_ch=PASS_Filter(impact_ch,"${params.GET_PASS_VARIANTS_PATH}")
}

workflow sv_cons {
    take:
        bam_ch
        cosmic_vcf_ch
    main:
        manta_fldr=file("${params.outdir}/Manta")

        if ( ! manta_fldr.isDirectory()) {
            result = manta_fldr.mkdir()
            println result ? "${manta_fldr}: OK" : "Cannot create directory: ${manta_fldr}"
        }

        manta_ch=Manta(bam_ch, "${params.ReferenceGenome}", "${params.ReferenceGenome}.fai")

        lumpy_fldr=file("${params.outdir}/Lumpy")

        if ( ! lumpy_fldr.isDirectory()) {
            result = lumpy_fldr.mkdir()
            println result ? "${lumpy_fldr}: OK" : "Cannot create directory: ${lumpy_fldr}"
        }

        lumpy_ch=Lumpy(bam_ch, "${params.ReferenceGenome}", "${params.ReferenceGenome}.fai", "${params.ex_bed_smoove}")

        cnvnator_fldr=file("${params.outdir}/CNVnator")

        if ( ! cnvnator_fldr.isDirectory()) {
            result = cnvnator_fldr.mkdir()
            println result ? "${cnvnator_fldr}: OK" : "Cannot create directory: ${cnvnator_fldr}"
        }

        cnvnator_ch=Cnvnator(bam_ch, "${params.ReferenceGenome}", "${params.ReferenceGenome}.fai")

        breakdancer_fldr=file("${params.outdir}/Breakdancer")

        if ( ! breakdancer_fldr.isDirectory()) {
            result = breakdancer_fldr.mkdir()
            println result ? "${breakdancer_fldr}: OK" : "Cannot create directory: ${breakdancer_fldr}"
        }

        breakdancer_ch=Breakdancer(bam_ch,"${params.ReferenceGenome}","${params.ReferenceGenome}.fai")

        for_mod_ch=manta_ch.join(lumpy_ch).join(cnvnator_ch)

        format_ch=VCF_format_mod(for_mod_ch)

        preprocess_ch=format_ch[1].join(format_ch[0]).join(format_ch[2]).join(breakdancer_ch)

        sv_type_ch=Consensus_preprocess(preprocess_ch, "${params.py_script_svtype}")

        ins_ch=Manta_ins_vcf(format_ch[1])

        survivor_del_dup_ch=sv_type_ch[0]
        survivor_inv_tra_ch=sv_type_ch[1]

        post_survivor_del_dup_ch=Survivor_del_dup(survivor_del_dup_ch)
        post_survivor_inv_tra_ch=Survivor_inv_tra(survivor_inv_tra_ch)

        fix_survivor_ch=post_survivor_del_dup_ch.join(post_survivor_inv_tra_ch)

        post_fix_ch=Fix_survivor_call(fix_survivor_ch)

        sv_fldr=file("${params.outdir}/SV")

        if ( ! sv_fldr.isDirectory()) {
            result = sv_fldr.mkdir()
            println result ? "${sv_fldr}: OK" : "Cannot create directory: ${sv_fldr}"
        }
	
        ind_ch=cosmic_vcf_ch.branch{
                gz:     it.size()==3
                no_ind: it.size()==2
        }

        vcf_germ_ch=Vcf_germ_ind(ind_ch.no_ind)

        ind_gz_ch=ind_ch.gz

        indexed_ch=ind_gz_ch.mix(vcf_germ_ch).ifEmpty(vcf_germ_ch.mix(ind_gz_ch))

        dup_ch=bam_ch.join(post_fix_ch[1]).join(indexed_ch)

        duphold_ch=Duphold_no_tra(dup_ch, "${params.ReferenceGenome}", "${params.ReferenceGenome}.fai") 

        postprocess_ch=duphold_ch.join(post_fix_ch[0]).join(ins_ch)

        post_ch=Postprocess_duphold(postprocess_ch)
    emit:
        post_ch
}   

workflow sv_annotation {
    take: postdup_vcf
    main:
        dir_SVAFotate=file("${params.outdir}/${params.SVAFOTATE_FOLDER}")
        if ( ! dir_SVAFotate.isDirectory() ) {
            result = dir_SVAFotate.mkdir()
            println result ? "${dir_SVAFotate}: OK" : "Cannot create directory: $dir_SVAFotate"
        }

        svafotate_ch=SVAFotate(postdup_vcf, file("${params.BED_SVAFotate}"))

        dir_AnnotSV=file("${params.outdir}/${params.ANNOTSV_FOLDER}")
        if ( ! dir_AnnotSV.isDirectory() ) {
            result = dir_AnnotSV.mkdir()
            println result ? "${dir_AnnotSV}: OK" : "Cannot create directory: $dir_AnnotSV"
        }

        annotsv_ch=AnnotSV(svafotate_ch)
    emit:
        annotsv_ch[0]
}

workflow sv_filt {
    take: 
        annoted_ch
        bam_ch
        melt_ch
    main:
        dir_tsvFiltering=file("${params.outdir}/${params.TSVFILTERING_FOLDER}")
        if ( ! dir_tsvFiltering.isDirectory() ) {
            result = dir_tsvFiltering.mkdir()
            println result ? "${dir_tsvFiltering}: OK" : "Cannot create directory: $dir_tsvFiltering"
        }

        filtering_ch=annoted_ch.join(bam_ch).join(melt_ch)

        TSV_filtering(filtering_ch, file("${params.BED_for_filtering}"), file("${params.TSV_DOSAGE}"), file("${params.outdir}/${params.TSVFILTERING_FOLDER}"))
}

workflow sv_filt_no_melt {
    take: 
        annoted_ch
        bam_ch
    main:
        dir_tsvFiltering=file("${params.outdir}/${params.TSVFILTERING_FOLDER}")
        if ( ! dir_tsvFiltering.isDirectory() ) {
            result = dir_tsvFiltering.mkdir()
            println result ? "${dir_tsvFiltering}: OK" : "Cannot create directory: $dir_tsvFiltering"
        }

        filtering_ch=annoted_ch.join(bam_ch)

        TSV_filtering_no_MELT(filtering_ch, file("${params.BED_for_filtering}"), file("${params.TSV_DOSAGE}"), file("${params.outdir}/${params.TSVFILTERING_FOLDER}"))
}

workflow melt{
    take:
        bam_ch
    main:
        dir_Melt=file("${params.outdir}/Melt")
        if ( ! dir_Melt.isDirectory() ) {
            result = dir_Melt.mkdir()
            println result ? "${dir_Melt}: OK" : "Cannot create directory: $dir_Melt"
        }
        
        melt_ch=Melt(bam_ch,
            file("/hpcshare/genomics/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta"),
            file("/hpcshare/genomics/references/gatk_bundle/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai"),
            "/opt/melt/add_bed_files/Hg38/Hg38.genes.bed",
            file("/hpcshare/genomics/bioinfo_parabricks/melt_mei_list.txt"))
        
        svafotate_ch=SVAFotate_Melt(melt_ch,
            file("${params.BED_SVAFotate}"))
        
        annotsv_ch=AnnotSV_Melt(svafotate_ch)

        concat_ch=Melt_concatenate(annotsv_ch)
    emit:
        concat_ch
}

workflow PB_germ {
    pb=pb_germline(fastq_ch)
    coverage(pb[0])
    add_pass(pb[1])
}

workflow Coverage_workflow {
    coverage(bam_ch)
}

workflow PB_germ_gvcf {
    pb=pb_germline_gvcf(fastq_ch)
    coverage(pb[0])
    add_pass(pb[1])
}

workflow AddPassHardFilter {
    add_pass(germ_snv_ch)
}

workflow CNVkit_workflow {
    cnvkit(bam_ch)
}

workflow SNV_annot {
    germline_snv_annotation(germ_snv_ch)
}


workflow SNV_AddPass_annot_filt{
    add_pass_ch=add_pass(germ_snv_ch)
    annotation_results=germline_snv_annotation(add_pass_ch)
    snv_filt(annotation_results)
}

workflow SNV_annot_filt {
    germline_snv_annotation(germ_snv_ch)
    snv_filt(germline_snv_annotation.out)   
}

workflow SNV_filt{
    snv_filt(annot_ch)
}

workflow SNV_filt_no_ONC{
    snv_filt_no_onc(annot_ch)
}

workflow SV_consensus {
    sv_ch=sv_cons(bam_ch,germ_snv_ch)
    annot_ch=sv_annotation(sv_ch)
    melt_ch=melt(bam_ch)
    sv_filt(annot_ch,bam_ch,melt_ch)
}

workflow SV_consensus_no_MELT {
    sv_ch=sv_cons(bam_ch,germ_snv_ch)
    annot_ch=sv_annotation(sv_ch)
    sv_filt_no_melt(annot_ch,bam_ch)
}

workflow SV_annot_filt {
    annot_ch=sv_annotation(germ_sv_ch)
    melt_ch=melt(bam_ch)
    sv_filt(annot_ch,bam_ch,melt_ch)
}

workflow SV_annot_filt_no_Melt {
    annot_ch=sv_annotation(germ_sv_ch)
    sv_filt_no_melt(annot_ch,bam_ch)
}

workflow Melt_workflow {
    melt(bam_ch)
}

workflow Pipeline_Germline_SNV_CNVkit{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    cnvkit(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    annotation_results=germline_snv_annotation(add_pass_ch)
    snv_filt(annotation_results)
}

workflow Pipeline_Germline_SNV{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    annotation_results=germline_snv_annotation(add_pass_ch)
    snv_filt(annotation_results)
}

workflow Pipeline_Germline_SNV_nofilt{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    annotation_results=germline_snv_annotation(add_pass_ch)
}

workflow Pipeline_Germline_CNVkit_SV{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    cnvkit(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    sv_ch=sv_cons(pb_results[0],add_pass_ch)
    annot_ch=sv_annotation(sv_ch)
    melt_ch=melt(pb_results[0])
    sv_filt(annot_ch,pb_results[0],melt_ch)
}

workflow Pipeline_Germline_SV{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    sv_ch=sv_cons(pb_results[0],add_pass_ch)
    annot_ch=sv_annotation(sv_ch)
    melt_ch=melt(pb_results[0])
    sv_filt(annot_ch,pb_results[0],melt_ch)
}

workflow Pipeline_Germline_SNV_nofilt_CNVkit_SV{
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    cnvkit(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    annotation_results=germline_snv_annotation(add_pass_ch)
    sv_ch=sv_cons(pb_results[0],add_pass_ch)
    annot_ch=sv_annotation(sv_ch)
    melt_ch=melt(pb_results[0])
    sv_filt(annot_ch,pb_results[0],melt_ch)
}

workflow {
    pb_results=pb_germline(fastq_ch)
    coverage(pb_results[0])
    cnvkit(pb_results[0])
    add_pass_ch=add_pass(pb_results[1])
    annotation_results=germline_snv_annotation(add_pass_ch)
    snv_filt(annotation_results)
    sv_ch=sv_cons(pb_results[0],add_pass_ch)
    annot_ch=sv_annotation(sv_ch)
    melt_ch=melt(pb_results[0])
    sv_filt(annot_ch,pb_results[0],melt_ch)
}

