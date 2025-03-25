process Impact_Filter {
input:
        tuple val(sample_id), path(vcf_impact_in)
output:
        tuple val(sample_id), path(vcf_impact_out)
script:
        sample_ID=vcf_impact_in.getSimpleName()
        vcf_impact_out="${sample_ID}_impact_filtered.vcf"
        """
        java -jar /usr/local/bin/snpEff/SnpSift.jar filter \
		 \"(ANN[*].IMPACT = 'HIGH') | (ANN[*].IMPACT = 'MODERATE' )\" \
		${vcf_impact_in} > ${vcf_impact_out}
        """
}
