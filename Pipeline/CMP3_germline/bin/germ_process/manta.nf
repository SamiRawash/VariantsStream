process Manta {
input:
        tuple val(sample_id), path(bam_path), path(ind_path) 
        path gh
        path gh_ind
output:
        tuple val(sample_id), path(manta_vcf_output)
shell:
        manta_vcf_output="${sample_id}_mantaGerm_INVconverted.vcf"
        '''
        echo "manta" !{bam_path}
        
        /swf/manta-1.6.0.centos6_x86_64/bin/configManta.py \
        --bam !{bam_path} \
        --referenceFasta !{gh} \
        --runDir ./

        ./runWorkflow.py -j !{task.cpus}

        /swf/manta-1.6.0.centos6_x86_64/libexec/convertInversion.py \
        /swf/manta-1.6.0.centos6_x86_64/libexec/samtools \
        !{gh} \
        ./results/variants/diploidSV.vcf.gz > \
        ./!{manta_vcf_output}
        '''
}