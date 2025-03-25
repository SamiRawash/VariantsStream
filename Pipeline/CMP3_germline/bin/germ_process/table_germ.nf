process Vcf2Table {
        tag "Vcf2Table"
input:
        tuple val(sample_id), path(final_vcf)
        path genelist_cantrans
output:
        path "${sample_ID}.tsv"
script:
        sample_ID=final_vcf.getSimpleName()
        """
        table -i ${final_vcf} -o ${final_vcf.getSimpleName()}.tsv -gl ${genelist_cantrans} -only -s
        """
}