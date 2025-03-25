process SnpEff_Annotation {
input:
        tuple val(sample_id), path(snpeff_in)
        val reference_genome
output:
        tuple val(sample_id), path(snpeff_out)
        path(report_snp)

script:
        sample_ID=snpeff_in.getSimpleName()
        snpeff_out="${sample_ID}_SnpEff.vcf"
        statis_snp="${sample_ID}_SnpEff.html"
        report_snp="${sample_ID}_SnpEff.csv"
        mem=task.memory.toString().split(" ")[0]
        """
        java -Xmx${mem}g -jar  /usr/local/bin/snpEff/snpEff.jar eff \
		-s "${statis_snp}" \
		-nodownload \
		-dataDir /home/dbs/ \
		-csvStats "${report_snp}" \
		${reference_genome} \
		${snpeff_in} > "${snpeff_out}"
        """
}

/* ***********************************************
******** VCF ANNOTATION - MultiQC SNPEFF  ********
************************************************ */

process MultiQC_SnpEff {
input :
        path(reports)
output :
        path "multiqc_report.html"
script:
        """
        multiqc .
        """
}
