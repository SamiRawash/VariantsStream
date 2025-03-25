process Breakdancer {
input:
        tuple val(sample_id), path(bam_path), path(ind_path)
        path gh
        path gh_ind
output:
        tuple val(sample_id), path(breakdancer_out)
shell:
        breakdancer_out="${sample_id}_S4.vcf"
        mem=task.memory.toString().split(" ")[0]
        '''
        echo "breakdancer" !{bam_path}
        /tools/SVE/bin/sve call !{bam_path} -r !{gh} -g hg38 -o ./ -t !{task.cpus} -M !{mem} -a breakdancer
        '''
}