process CosmicAnnotation {
input:
        tuple val(sample_id), path(input)
        path database
        path database_index
        val reference_genome
output:
        tuple val(sample_id), path(vcf_COSMIC_out)
script:
        sample_ID=input.getSimpleName()
        vcf_COSMIC_out="${sample_ID}_COSMIC.vcf"
        mem=task.memory.toString().split(" ")[0]
        """
        java -Xmx"${mem}"g -jar /usr/local/bin/snpEff/SnpSift.jar \
		annotate \
		"$database" \
		"$input"  > "$vcf_COSMIC_out"
        """
}
