process Lumpy {
input:
        tuple val(sample_id), path(bam_path), path(ind_path)
        path gh
        path gh_ind
        path ex
output:
        tuple val(sample_id), path(lumpy_out)
script:
        lumpy_out="${sample_id}-smoove.genotyped.vcf.gz"
        """
        echo "lumpy" ${bam_path}
        smoove call \
	--fasta ${gh} \
	--name ${sample_id} \
	-x \
	--genotype \
	--exclude ${ex} \
	-o ./ \
	${bam_path}
        """
}