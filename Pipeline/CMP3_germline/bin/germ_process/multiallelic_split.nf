process MultiAllelicSplit {
input:
        tuple val(sample_id), path(vcf_filtered), file(filtered_ind)
output:
        tuple val(sample_id), path(multiallelicsplitted)
script:
        sample_ID=filtered_ind.getSimpleName()
        multiallelicsplitted="${sample_ID}_MultiallelicSplit.vcf"
        """
        bcftools norm -m - ${vcf_filtered} > ${multiallelicsplitted}
        """
}