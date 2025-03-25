process AnnotSV {
input:
        tuple val(sample_id), path(vcf_dup_final)
output:
        tuple val(sample_id), path("./${sample_ID}_consensusSV_annotsv.tsv")
        tuple val(sample_id), path("./${sample_ID}_consensusSV_annotsv.unannotated.tsv")
script:
        sample_ID=vcf_dup_final.getSimpleName()
        """
        echo ${vcf_dup_final}

        /swf/AnnotSV/bin/AnnotSV -SVinputFile ${vcf_dup_final} \
                -outputFile ${sample_ID}_consensusSV_annotsv \
                -bedtools /swf/bedtools2/bin/bedtools \
                -bcftools /swf/bin/bcftools \
                -genomeBuild GRCh38
        ls

        mv ./*_AnnotSV/*.tsv ./
        """
}
