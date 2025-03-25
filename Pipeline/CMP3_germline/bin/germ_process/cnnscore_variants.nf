process CNNScoreVariants {
	executor 'pbspro'
        queue "gpunodes"
        memory "150 GB"
        cpus 16
        time "5h"
        tag 'CNNScoreVariants'
        publishDir "${baseDir}/${idSample}" 
input:
        tuple val(sample_id), file(bam), file(bai) 
        tuple val(sample_id), file(vcf) 
        path input_ref 
        path ind_ref
        val pb_exe
output:
        tuple val(sample_id), file(vcfcnn)
script:
        splitted=sample_id.split("_")
	idSample=splitted[2]
        vcfcnn="${sample_id}_CNN.vcf"
        
        """
        time "${pb_exe}" cnnscorevariants --ref "${input_ref}" --in-bam "${bam}" --in-vcf "${vcf}" --out-vcf "${vcfcnn}" --num-gpus 2
        """
}