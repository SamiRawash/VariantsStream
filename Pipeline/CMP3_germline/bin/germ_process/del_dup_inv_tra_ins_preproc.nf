process Manta_ins_vcf {
input:
        tuple val(sample_id), path(manta_vcf)
output:
        tuple val(sample_id), path(manta_ins_out), path("${manta_ins_out}.tbi")
script:
        manta_ins_out="${sample_id}_INS.vcf.gz"
        """                                                 
        echo "INS" ${manta_vcf} 

        echo ${sample_id} > ${sample_id}_NAME_sample.txt
        
        #Creo un VCF con solo le INS di MANTA

        grep -E '^#|SVTYPE=INS' ${manta_vcf} > ${sample_id}_mantaOK_toolName_onlyINS.vcf
        bcftools reheader -s  ${sample_id}_NAME_sample.txt ${sample_id}_mantaOK_toolName_onlyINS.vcf > ${sample_id}_INS.vcf
        bgzip -f -c ${sample_id}_INS.vcf > ${manta_ins_out}
        tabix ${manta_ins_out}
        """
}

process Consensus_preprocess {
input:
        tuple val(sample_id), path(manta_vcf), path(lumpy_vcf), path(cnvnator_vcf), path(breakdancer_vcf)
        val python_script
output:
        tuple val(sample_id), path("${manta_out}_DEL.vcf"), path("${manta_out}_DUP.vcf"),\
                path("${lumpy_out}_DEL.vcf"), path("${lumpy_out}_DUP.vcf"),\
                path("${cnvnator_out}_DEL.vcf"), path("${cnvnator_out}_DUP.vcf"),\
                path("${breakdancer_out}_DEL.vcf"), path("${breakdancer_out}_DUP.vcf")
        tuple val(sample_id), path("${manta_out}_INV.vcf"), path("${manta_out}_TRA.vcf"),\
                path("${lumpy_out}_INV.vcf"), path("${lumpy_out}_TRA.vcf"),\
                path("${breakdancer_out}_INV.vcf"), path("${breakdancer_out}_TRA.vcf")

script:
        manta_out="${sample_id}_mantaOK_toolName"
        lumpy_out="${sample_id}_S18OK_toolName"
        cnvnator_out="${sample_id}_S1O_toolName"
        breakdancer_out="${sample_id}_S4"
        """
        python3 ${python_script} --input ${manta_vcf} --sample_id ${manta_out}

        python3 ${python_script} --input ${lumpy_vcf} --sample_id ${lumpy_out}

        python3 ${python_script} --input ${cnvnator_vcf} --sample_id ${cnvnator_out}

        python3 ${python_script} --input ${breakdancer_vcf} --sample_id ${breakdancer_out}
        """
}