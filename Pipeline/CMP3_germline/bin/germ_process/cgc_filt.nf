process CGC_Filter {
input:
        tuple val(sample_id), path(vcf_CGC_in)
        path cgc_script
        path cgc_file
output:
        tuple val(sample_id), path(vcf_CGC_out)
script:
	sample_ID=vcf_CGC_in.getSimpleName()
        vcf_CGC_out="${sample_ID}_CGC.vcf"
        """
        #/home/analysis/scripts/cancer_gene_census_filter_python_script/cancer_gene_census_filter.py
        python3 "${cgc_script}" \
                -g ${cgc_file} \
                -v ${vcf_CGC_in} > ${vcf_CGC_out}
	"""
}
