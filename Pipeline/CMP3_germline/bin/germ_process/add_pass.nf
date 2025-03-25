process AddPass_VariantFiltration {
        publishDir "${outdir}/${sample}"
//        publishDir "${outdir}/"
input:
        tuple val(sample_id), file(vcf)
        val qd_value //2.0
	val qual_value //30.0
	val sor_value //3.0
	val fs_value //60.0
	val mq_value //40.0
	val mqranksum_value //-12.5
	val readposranksum_value //-8.0
        val outdir
output:
        tuple val(sample_id), path(vcf_out), path("${vcf_out}.tbi")

script:
        //template "VCF_filtration_DSL1.sh"
        splitted=sample_id.split("_")
        sample=splitted[2]
        vcf_out="${sample_id}_PassHF.vcf.gz"
	"""
	gatk VariantFiltration \
    		-V ${vcf} \
		--filter-expression "QD < ${qd_value}" --filter-name "QD2" \
		--filter-expression "QUAL < ${qual_value}" --filter-name "QUAL30" \
		--filter-expression "SOR > ${sor_value}" --filter-name "SOR3" \
		--filter-expression "FS > ${fs_value}" --filter-name "FS60" \
		--filter-expression "MQ < ${mq_value}" --filter-name "MQ40" \
		--filter-expression "MQRankSum < ${mqranksum_value}" --filter-name "MQRankSum-12.5" \
		--filter-expression "ReadPosRankSum < ${readposranksum_value}" --filter-name "ReadPosRankSum-8" \
		-O ${vcf_out}
	"""
}

process AddPass_VariantFiltration_gz {
        publishDir "${outdir}/${sample}"
//		publishDir "${outdir}/"
input:
        tuple val(sample_id), file(vcf), file(ind)
        val qd_value //2.0
	val qual_value //30.0
	val sor_value //3.0
	val fs_value //60.0
	val mq_value //40.0
	val mqranksum_value //-12.5
	val readposranksum_value //-8.0
        val outdir
output:
        tuple val(sample_id), path(vcf_out), path("${vcf_out}.tbi")

script:
        //template "VCF_filtration_DSL1.sh"
        splitted=sample_id.split("_")
        sample=splitted[2]
        vcf_out="${sample_id}_PassHF.vcf.gz"
	"""
	gatk VariantFiltration \
    		-V ${vcf} \
		--filter-expression "QD < ${qd_value}" --filter-name "QD2" \
		--filter-expression "QUAL < ${qual_value}" --filter-name "QUAL30" \
		--filter-expression "SOR > ${sor_value}" --filter-name "SOR3" \
		--filter-expression "FS > ${fs_value}" --filter-name "FS60" \
		--filter-expression "MQ < ${mq_value}" --filter-name "MQ40" \
		--filter-expression "MQRankSum < ${mqranksum_value}" --filter-name "MQRankSum-12.5" \
		--filter-expression "ReadPosRankSum < ${readposranksum_value}" --filter-name "ReadPosRankSum-8" \
		-O ${vcf_out}
	"""
}







process AddPass_FilterVariantTranches {
        container "${gatk_cont}"
input:
        tuple val(sample_id),file(input)
        path REF_HAPMAP_FILTERGERMLINE 
        path REF_HAPMAP_FILTERGERMLINE_index
        path REF_MILLS_FILTERGERMLINE
        path REF_MILLS_FILTERGERMLINE_index
        val gatk_cont
output:
        tuple val(sample_id), path(output), path("${output}.tbi")
script:
        sample_ID=input.getSimpleName()
        output="${sample_ID}_PassVT.vcf.gz"
        """
        gatk FilterVariantTranches -V ${input} --resource "${REF_HAPMAP_FILTERGERMLINE}" --resource "${REF_MILLS_FILTERGERMLINE}" --info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4 --invalidate-previous-filters -O ${output}
        """
}
