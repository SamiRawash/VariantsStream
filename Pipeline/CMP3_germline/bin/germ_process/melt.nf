process Melt {
input:
    tuple val(sample_id), path(BAM), path(BAI)
    file GENhg38
    file GENhg38_ind
    val hg38_genes_bed
    file mei_list
output:
    tuple val(sample_id), path("ALU.final_comp.vcf"), path("LINE1.final_comp.vcf"), path("HERVK.final_comp.vcf"), path("SVA.final_comp.vcf")  
script:
    """
    java -jar /opt/melt/MELT.jar Preprocess -bamfile ${BAM} -h ${GENhg38}
    
    ls 
    
    java -jar /opt/melt/MELT.jar Single -k -h ${GENhg38} \
        -bamfile ${BAM} \
        -n ${hg38_genes_bed} \
        -t ${mei_list} \
        -w ./
    
    ls
    """
}

process SVAFotate_Melt {
input:
    tuple val(sample_id), path(ALU), path(LINE1), path(HERVK), path(SVA)
    path bedSVAFotate 
output:
    tuple val(sample_id), path("ALU_MeltSVAFotate.vcf"), path("LINE1_MeltSVAFotate.vcf"), path("HERVK_MeltSVAFotate.vcf"), path("SVA_MeltSVAFotate.vcf")
shell:
    '''
    for mei in ALU LINE1 HERVK SVA; do
        VCF=${mei}.final_comp.vcf
        vcf_annotated=./${mei}_MeltSVAFotate.vcf
        svafotate annotate -v ${VCF} \
            -o ${vcf_annotated} \
            -b !{bedSVAFotate} \
            -a best -f 0.5
    done
    '''
}

process AnnotSV_Melt {
input:
    tuple val(sample_id), path(ALU), path(LINE1), path(HERVK), path(SVA)
output:
    tuple val(sample_id), path("ALU_MeltSVAFotate_melt_annotsv.tsv"), path("LINE1_MeltSVAFotate_melt_annotsv.tsv"), path("HERVK_MeltSVAFotate_melt_annotsv.tsv"), path("SVA_MeltSVAFotate_melt_annotsv.tsv")
shell:
    '''
    for mei in ALU LINE1 HERVK SVA; do
        VCF=${mei}_MeltSVAFotate.vcf
        vcf_annotated=./${mei}_MeltSVAFotate_melt_annotsv.tsv
        /swf/AnnotSV/bin/AnnotSV -SVinputFile ${VCF} \
        -outputFile ${vcf_annotated} \
        -bedtools /swf/bedtools2/bin/bedtools \
        -bcftools /swf/bin/bcftools \
        -genomeBuild GRCh38
    done
    '''
}

process Melt_concatenate{
input:
    tuple val(sample_id), path(ALU), path(LINE1), path(HERVK), path(SVA)
output:
    tuple val(sample_id), path("${sample_id}_MeltSVAFotate_melt_annotsv.tsv")
script:
    """
    concatenate -i ${ALU},${LINE1},${HERVK},${SVA} \
        -o ${sample_id}_MeltSVAFotate_melt_annotsv.tsv
    """
}