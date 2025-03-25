process Cnvnator {
input:
        tuple val(sample_id), path(bam_path), path(ind_path)
        path gh
        path gh_ind
output:
        tuple val(sample_id), path(cnvnator_out)
shell:
        cnvnator_out="${sample_id}_S10.vcf"
        mem=task.memory.toString().split(" ")[0]
        '''
        echo "cnvnator" !{bam_path}
        /tools/SVE/bin/sve call !{bam_path} -r !{gh} -g hg38 -o ./ -t !{task.cpus} -M !{mem} -a cnvnator
        '''
}
