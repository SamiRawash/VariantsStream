process AF_Filter {
input:
        tuple val(sample_id), path(vcf_AF_in)
       	val af
	val af_value
output:
        tuple val(sample_id), path(vcf_AF_out)
script:
        sample_ID=vcf_AF_in.getSimpleName()
        vcf_AF_out="${sample_ID}_AF_filtered.vcf"
        """
        java -jar /usr/local/bin/snpEff/SnpSift.jar \
		filter  \"(${af} <= ${af_value}) | !(exists ${af})\" \
		${vcf_AF_in} > ${vcf_AF_out}
        """
}
