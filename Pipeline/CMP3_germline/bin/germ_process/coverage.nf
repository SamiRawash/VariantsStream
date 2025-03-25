process Coverage {
input:
    tuple val(sample_id), file(bam_file),file(bai_file)
output:
	file coverage_out 

script:
    coverage_out="coverage_${sample_id}.txt"
	"""
    samtools coverage "${bam_file}" --output "${coverage_out}"  
    """
}
