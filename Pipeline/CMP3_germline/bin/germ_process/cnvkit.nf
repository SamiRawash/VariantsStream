process CNVkit {
	publishDir "${outdir}/${sample}" 
//	publishDir "${outdir}/"
input:
	tuple val(sample_id), file(bam_file), file(bai_file) 
	path reference 
	val CNVKIT_EXEC
	val outdir
output:
	tuple val(sample_id),file("${sample_id}.cns") ,file("${sample_id}.cnr"),file("${sample_id}.call.cns")
script:
	splitted=sample_id.split("_")
    	sample=splitted[2]
	"""
	${CNVKIT_EXEC} batch "${bam_file}" -r "${reference}" --scatter --diagram -p ${params.CNVKIT_CPUS}
	"""
}

process CNVkit_filter {
	publishDir "${outdir}/${sample}"
//	publishDir "${outdir}/"
input:
	tuple val(sample_id), file(file_cns),file(file_cnr),file(file_call_cns)
	path gene_to_filter 
	val CNVKIT_FILTERING_SCRIPT
	val CNVKIT_FILTERING_MODIFY_SCRIPT
	val outdir
output:
	tuple val(sample_id), file(file_cns),file(file_cnr),file(file_call_cns),file("${sample_id}_filtered_log2_genes.txt"),file("${sample_id}_modified.cnr"),file("${sample_id}_filtered_log2_genes_call.txt")
script:
	splitted=sample_id.split("_")
    	sample=splitted[2]
	"""
	python3 "${CNVKIT_FILTERING_SCRIPT}" -c "${file_call_cns}" -g "${gene_to_filter}" -o "${sample_id}_filtered_log2_genes_call.txt"
	python3 "${CNVKIT_FILTERING_SCRIPT}" -c "${file_cns}" -g "${gene_to_filter}" -o "${sample_id}_filtered_log2_genes.txt"
	python3 "${CNVKIT_FILTERING_MODIFY_SCRIPT}" -f "${file_cnr}" -o "${sample_id}_modified.cnr"
	"""
}

process CNVkit_plot {
	publishDir "${outdir}/${sample}", mode: 'move', overwrite: true
//	publishDir "${outdir}/${sample}", mode: 'move', overwrite: true
input:
	tuple val(sample_id), file(cns),file(cnr),file(altereted_genes),file(all_genes),file(log2_call)
	path rscript_plot 
	path centromere_position
	val outdir
output:
	path "${sample_id}_chr*.png"
script:
	splitted=sample_id.split("_")
    	sample=splitted[2]
	"""
	Rscript --vanilla "${rscript_plot}" "${cnr}" "${cns}" "${all_genes}" "${altereted_genes}" "${centromere_position}"
	"""
}
