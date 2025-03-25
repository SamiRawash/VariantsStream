process VCF_format_mod {
input:
        tuple val(sample_id), path(manta_vcf), path(lumpy_vcf), path(cnvnator_vcf)
output:
        tuple val(sample_id), path(lumpy_form_out)
        tuple val(sample_id), path(manta_form_out)
        tuple val(sample_id), path(cnvnator_form_out)
script:
        lumpy_form_out="${sample_id}_S18OK_toolName.vcf"
        manta_form_out="${sample_id}_mantaOK_toolName.vcf"
        cnvnator_form_out="${sample_id}_S10_toolName.vcf"
        """
        echo "mod" ${manta_vcf} ${lumpy_vcf} ${cnvnator_vcf}

        sed 's/SVTYPE=BND/SVTYPE=TRA/' ${manta_vcf} > ${sample_id}_mantaOK.vcf 
        zcat ${lumpy_vcf} | sed 's/SVTYPE=BND/SVTYPE=TRA/'  > ${sample_id}_S18OK.vcf

        ####Creo dei files con il nome per il campo format da andare a modificare successivamente nel VCF.
        echo "${sample_id}_LUMPY" > ${sample_id}_lumpy.txt
        echo "${sample_id}_MANTA" > ${sample_id}_manta.txt
        echo "${sample_id}_CNVNATOR" > ${sample_id}_cnvnator.txt
        echo "${sample_id}_BREAKDANCER" > ${sample_id}_breakdancer.txt

        ####Modifico il nome del campo format nel VCF - serve successivamente quando andr� ad usare lo script fix_SURVIVORgenotypes.R
        #Per breakdancer non lo faccio perchè non ha campo format
        echo "reheader lumpy"
        bcftools reheader -s ${sample_id}_lumpy.txt ${sample_id}_S18OK.vcf > ${lumpy_form_out}
        echo "reheader manta"
        bcftools reheader -s ${sample_id}_manta.txt ${sample_id}_mantaOK.vcf > ${manta_form_out}
        echo "reheader cnvnator"
        bcftools reheader -s ${sample_id}_cnvnator.txt ${cnvnator_vcf} > ${cnvnator_form_out}
        """
}