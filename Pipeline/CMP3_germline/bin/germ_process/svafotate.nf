process SVAFotate {
input:
    tuple val(sample_id), path(vcf_dup_final) //from vcfs_consensus 
    path bedSVAFotate //from params.BED_SVAFotate
output:
    tuple val(sample_id), path(SVAFotate_out) //into SVAFotate
script:
    SVAFotate_out="${vcf_dup_final.getSimpleName()}_SVAFotate.vcf"
    """
    svafotate annotate -v ${vcf_dup_final} -o ${SVAFotate_out} -b ${bedSVAFotate} -a best -f 0.5
    """
}
