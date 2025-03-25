process PipelineGermline { 
     tag 'ParabricksGemlineLine'
     //publishDir "${outdir}/"
     publishDir "${outdir}/${sample}"
     containerOptions "--nv"

input:
    tuple val(sample_id),file(fq1),file(fq2) 
	path input_ref 
	path input_knownsites 
	val outdir
output:
	tuple val(sample_ID), file(bam), file("${bam}.bai") 
    tuple val(sample_ID), file(vcf) 
    file rcl 
    file dpl 
shell:
	splitted=sample_id.split("_")
	comp=splitted[0]
    project=splitted[1]
    sample=splitted[2]
    sample_ID="${comp}_${project}_${sample}"
    pl="${params.sequencing_platform}"
	lb="${params.sequencing_lib}"
	bam="${comp}_${project}_${sample}.bam"
	vcf="${comp}_${project}_${sample}.vcf"
	rcl="${comp}_${project}_${sample}_Recal.txt"
	dpl="${comp}_${project}_${sample}_DuplMetr.txt"
	'''
	PU=$( zless !{fq1} | head -1 | awk -v sample=!{sample} '{split($1,vect1,":") ; split($2,vect2,":") ; printf("%s:%s:%s.%s.%s", substr(vect1[1],2), vect1[2],vect1[3],vect1[4],sample) }' )
	LB=$( zless !{fq1} | head -1 | awk -v lb=!{lb} '{split($2,vect2,":") ; printf("%s:INDX_%s", lb, vect2[4]) }' )
	time pbrun germline --ref "!{input_ref}" --in-fq "!{fq1}" "!{fq2}" --read-group-sm "!{sample}" --read-group-lb "${LB}" --read-group-pl "!{pl}" --read-group-id-prefix "${PU}" --knownSites "!{input_knownsites}" --out-bam "!{bam}" --out-variants "!{vcf}" --out-recal-file "!{rcl}" --out-duplicate-metrics "!{dpl}" --num-gpus !{params.pb_ngpus}
	'''
}

process PipelineGermline_gvcf {
        tag 'ParabricksGemlineLine'
//    publishDir "${outdir}/"
    publishDir "${outdir}/${sample}"
    containerOptions "--nv"
input:
    tuple val(sample_id),file(fq1),file(fq2)
        path input_ref
        path input_knownsites
        val outdir
output:
        tuple val(sample_ID), file(bam), file("${bam}.bai")
    tuple val(sample_ID), file(vcf), file("${vcf}.tbi")
    file rcl
    file dpl
shell:
        splitted=sample_id.split("_")
        comp=splitted[0]
    project=splitted[1]
    sample=splitted[2]
    sample_ID="${comp}_${project}_${sample}"
    pl="${params.sequencing_platform}"
        lb="${params.sequencing_lib}"
        bam="${comp}_${project}_${sample}.bam"
        vcf="${comp}_${project}_${sample}.g.vcf.gz"
        rcl="${comp}_${project}_${sample}_Recal.txt"
        dpl="${comp}_${project}_${sample}_DuplMetr.txt"
        '''
        PU=$( zless !{fq1} | head -1 | awk -v sample=!{sample} '{split($1,vect1,":") ; split($2,vect2,":") ; printf("%s:%s:%s.%s.%s", substr(vect1[1],2), vect1[2],vect1[3],vect1[4],sample) }' )
        LB=$( zless !{fq1} | head -1 | awk -v lb=!{lb} '{split($2,vect2,":") ; printf("%s:INDX_%s", lb, vect2[4]) }' )
        time pbrun germline --ref "!{input_ref}" --in-fq "!{fq1}" "!{fq2}" \
                --read-group-sm "!{sample}" \
                --read-group-lb "${LB}" \
                --read-group-pl "!{pl}" \
                --read-group-id-prefix "${PU}" \
                --knownSites "!{input_knownsites}" \
                --out-bam "!{bam}" \
                --out-variants "!{vcf}" \
                --gvcf \
                --out-recal-file "!{rcl}" \
                --out-duplicate-metrics "!{dpl}" \
                --num-gpus !{params.pb_ngpus}
        '''
}

process Haplotype_gvcf {
    tag 'ParabricksGemlineLine'
    publishDir "${outdir}/${sample}"
    //publishDir "${outdir}/"    
    containerOptions "--nv"
input:
    tuple val(sample_id), path(bam_path), path(ind_path), 
    path input_ref
    path input_knownsites
    val outdir
output:
    tuple val(sample_id), file(vcf), file("${vcf}.tbi")
shell:
    vcf="${sample_id}.g.vcf.gz"
    '''
    time pbrun haplotypecaller --ref "!{input_ref}" \
            --in-bam "!{bam}" \
            --out-variants "!{vcf}" \
            --in-recal-file "!{rcl}" \
            --gvcf \
            --num-gpus !{params.pb_ngpus}
    '''
}
