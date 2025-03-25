process IncidentalFindings_filter {
        tag "Filter4Report"
input:
        tuple val(sample_id), path(vcf)
output:
        tuple val(sample_id), path("${sample_ID}_PASS4Report.vcf")
script:
        sample_ID=vcf.getSimpleName()
        """
        incidental_findings -o ${sample_ID}_PASS4Report.vcf -i ${vcf} 
        """
}

process Vcf2Table_IF {
        tag "Vcf2Table_IF"
input:
        tuple val(sample_id),path(IF_vcf)
        path genelist_cantrans
output:
        path "${sample_ID}.tsv"
script:
        sample_ID=IF_vcf.getSimpleName()
        """
        table -i ${IF_vcf} -o ${IF_vcf.getSimpleName()}.tsv -gl ${genelist_cantrans} -only -s -canonical
        """
}
