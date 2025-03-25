process Vcf_germ_ind {
input:
        tuple val(sample_id), path(germ)
output:
        tuple val(sample_ID), path("${germ}.gz"), path("${germ}.gz.tbi")
script:
        splitted=sample_id.split("_")
        comp=splitted[0]
        project=splitted[1]
        sample=splitted[2]
        sample_ID="${comp}_${project}_${sample}"
        """
	file ${germ}
	htsfile ${germ}
	bgzip -f -c ${germ} > ${germ}.gz
	ls
	tabix ${germ}.gz
        """
}

process Duphold_no_tra {
input:
        tuple val(sample_id), path(bam_path), path(bam_ind), path(vcf_fixGenotype_fixCIPOS_noTRA), path(vcf_germ), path(ind_germ)
        path gh
        path gh_ind
output:
        tuple val(sample_id), path(duphold_out)
script:
        duphold_out="${sample_id}-smoove.genotyped_Duphold.vcf"
        """
        echo ${bam_path} ${vcf_fixGenotype_fixCIPOS_noTRA} ${vcf_germ}
        echo "genotyping"
        smoove genotype \
                --fasta ${gh} \
                --name ${sample_id} \
                -d -x \
                -o ./ \
                --vcf ${vcf_fixGenotype_fixCIPOS_noTRA} \
                ${bam_path}
        
        if [ -s ${vcf_germ} ]; then
                smoove duphold \
                        --fasta ${gh} \
                        --vcf ./${sample_id}-smoove.genotyped.vcf.gz \
                        --snps ${vcf_germ} \
                        -o ./${duphold_out}.gz \
                        ${bam_path}   
        else
                mv ./${sample_id}-smoove.genotyped.vcf.gz ./${sample_id}-smoove.genotyped_Duphold.vcf.gz
        fi
        gunzip ./${duphold_out}.gz
        """
}

process Postprocess_duphold {
input:
        tuple val(sample_id), path(duphold_vcf), path(vcf_TRA_fixed), path(vcf_TRA_fixed_ind), path(vcf_manta_ins), path(vcf_manta_ins_ind)
output:
        tuple val(sample_id), path(postdup_out)
script: 
        postdup_out="${sample_id}_consensusFiltered_final.vcf"
        """
        echo ${duphold_vcf} ${vcf_TRA_fixed} ${vcf_manta_ins}

        echo ${sample_id} > ${sample_id}_NAME_sample.txt
        bcftools reheader -s  ${sample_id}_NAME_sample.txt ${duphold_vcf} > ${sample_id}-smoove.genotyped_Duphold_nameOK.vcf
        #bcftools index -t ${sample_id}-smoove.genotyped_Duphold_nameOK.vcf
        echo -e "\n ########## RUNNING FILTERS ON DUPHOLD ##########"

        Rscript /opt/filterSVByFC.R \
                ${sample_id}-smoove.genotyped_Duphold_nameOK.vcf \
                ${sample_id}-smoove.genotyped_Duphold_nameOK_filtered.vcf.gz

        Rscript /opt/rename_SVIds.R \
                ${sample_id}-smoove.genotyped_Duphold_nameOK_filtered.vcf.gz \
                ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID.vcf.gz #OK

        echo -e "\n ########## RUNNING FILTERS ON INFO TAGS ##########"
        zcat ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID.vcf.gz \
                | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" \
                >${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID_BCFTOOLS.vcf

        ####Mergiare VCF con DEL, DUP, INS che ha fatto anche smoove genotype CON quello contenente solo le TRA
        bgzip -f -c ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID_BCFTOOLS.vcf > ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID_BCFTOOLS.vcf.gz
        bcftools index -t ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID_BCFTOOLS.vcf.gz

        echo -e "\n ########## RUNNING BCFTOOLS CONCAT to MERGE VCF with DEL DUP INV already filtered with DUPHOLD with VCF with TRA ##########"
        bcftools concat -a -O v -o ${sample_id}_consensusFiltered_final.vcf \
                ${sample_id}-smoove.genotyped_nameOK_dupholdFiltered_renameID_BCFTOOLS.vcf.gz \
                ${vcf_TRA_fixed} \
                ${vcf_manta_ins}
        """
}
