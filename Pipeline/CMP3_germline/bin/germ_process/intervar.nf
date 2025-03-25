process InterVar_Filter {
input:
        tuple val(sample_id), path(vcf_InterVar_in)
        val InterVar_py_script
output:
        tuple val(sample_id), path(vcf_INTERVAR_out)
script:
        sample_ID=vcf_InterVar_in.getSimpleName()
        vcf_INTERVAR_out="${sample_ID}_InterVar.vcf"
        """
        python3 ${InterVar_py_script} -v ${vcf_InterVar_in} -o ${vcf_INTERVAR_out}
        """
}