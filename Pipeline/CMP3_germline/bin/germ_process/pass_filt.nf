process PASS_Filter {
input:
        tuple val(sample_id), path(vcf_PASS_in)
        val get_pass_variants_path
output:
        tuple val(sample_id), path(vcf_PASS_out)
script:
	sample_ID=vcf_PASS_in.getSimpleName()
        vcf_PASS_out="${sample_ID}_PASS_filtered.vcf"
        """
	${get_pass_variants_path} ${vcf_PASS_in} > ${vcf_PASS_out}
        """
}
