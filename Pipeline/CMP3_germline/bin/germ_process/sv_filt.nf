process TSV_filtering {
input:
    tuple val(sample_id), file(tsv),file(bam_file),file(bai_file), file(melt_tsv)
    path filtering_bed 
    path dosage_tsv
    path output_dir
shell:
	splitted="${tsv.getSimpleName()}".split("_")
    sample=splitted[2]
    '''
	filterSV -i !{tsv} \
        -b !{bam_file} \
        -g !{filtering_bed} \
		-f \
        -d !{dosage_tsv} \
        -m !{melt_tsv} \
        -o !{output_dir}/!{sample}/ 
    
	ls
	'''
}

process TSV_filtering_no_MELT {
input:
    tuple val(sample_id), file(tsv),file(bam_file),file(bai_file)
    path filtering_bed 
    path dosage_tsv
    path output_dir
shell:
	splitted="${tsv.getSimpleName()}".split("_")
    sample=splitted[2]
    '''
	filterSV -i !{tsv} \
        -b !{bam_file} \
        -g !{filtering_bed} \
		-f \
        -d !{dosage_tsv} \
        -o !{output_dir}/!{sample}/  
    
	ls
	'''
}