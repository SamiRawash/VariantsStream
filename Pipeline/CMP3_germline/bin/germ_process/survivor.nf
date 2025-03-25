process Survivor_del_dup {
input:
        tuple val(sample_id), path(manta_vcf_dup), path(manta_vcf_del), path(lumpy_vcf_dup), path(lumpy_vcf_del), path(cnvnator_vcf_dup), path(cnvnator_vcf_del), path(breakdancer_vcf_dup), path(breakdancer_vcf_del)
output:
        tuple val(sample_id), path(survivor_del_out), path(survivor_dup_out)
shell:
        survivor_del_out="${sample_id}_DEL_consensus.vcf"
        survivor_dup_out="${sample_id}_DUP_consensus.vcf"
        '''
        echo !{manta_vcf_dup} !{lumpy_vcf_dup} !{cnvnator_vcf_dup} !{breakdancer_vcf_dup} !{manta_vcf_del} !{lumpy_vcf_del} !{cnvnator_vcf_del} !{breakdancer_vcf_del}

        breakpoint_dist=1000
        min_tool_calls=2
        use_type=1
        use_strand=1
        dist_based=0
        min_sv_size=0

        for type in DEL DUP; do
                #Creo file con lista dei VCFs contenenti le DEL di ogni tool.
                ls *_${type}.vcf > !{sample_id}_sample_file_${type}.txt
                echo "########## DOING CONSENSUS - SURVIVOR: $type"
                SURVIVOR merge \
                        !{sample_id}_sample_file_${type}.txt \
                        ${breakpoint_dist} \
                        ${min_tool_calls} \
                        ${use_type} \
                        ${use_strand} \
                        ${dist_based} \
                        ${min_sv_size} \
                                !{sample_id}_${type}_consensus.vcf
                number=$(grep -v "^#" !{sample_id}_${type}_consensus.vcf | wc -l)
                echo -e "SURVIVOR: Number of consensus in VCF obtained: ${number} \n"
        done
        '''
}

process Survivor_inv_tra {
input:
        tuple val(sample_id), path(manta_vcf_inv), path(manta_vcf_tra), path(lumpy_vcf_inv), path(lumpy_vcf_tra), path(breakdancer_vcf_inv), path(breakdancer_vcf_tra)
output:
        tuple val(sample_id), path(survivor_tra_out), path(survivor_inv_out)
shell:
        survivor_tra_out="${sample_id}_TRA_consensus.vcf"
        survivor_inv_out="${sample_id}_INV_consensus.vcf"
        '''
        echo !{manta_vcf_inv} !{lumpy_vcf_inv} !{breakdancer_vcf_inv} !{manta_vcf_tra} !{lumpy_vcf_tra} !{breakdancer_vcf_tra}

        breakpoint_dist=1000
        min_tool_calls=2
        use_type=1
        use_strand=1
        dist_based=0
        min_sv_size=0

        for type in TRA INV; do
                #Creo file con lista dei VCFs contenenti le DEL di ogni tool.
                ls *_${type}.vcf > !{sample_id}_sample_file_${type}.txt
                echo "########## DOING CONSENSUS - SURVIVOR: $type"
                SURVIVOR merge \
                        !{sample_id}_sample_file_${type}.txt \
                        ${breakpoint_dist} \
                        ${min_tool_calls} \
                        ${use_type} \
                        ${use_strand} \
                        ${dist_based} \
                        ${min_sv_size} \
                                !{sample_id}_${type}_consensus.vcf
                number=$(grep -v "^#" !{sample_id}_sample_file_${type}_consensus.vcf | wc -l)
                echo -e "SURVIVOR: Number of consensus in VCF obtained: ${number} \n"
        done
        '''
}

process Fix_survivor_call {
input:
        tuple val(sample_id), path(vcf_DEL_cons), path(vcf_DUP_cons), path(vcf_TRA_cons), path(vcf_INV_cons)
output:
        tuple val(sample_id), path(fix_tra_out), path("${fix_tra_out}.tbi")
        tuple val(sample_id), path(fix_merge_out)
shell:
        fix_tra_out="${sample_id}_TRA_consensus_sorted_fixed_cipos.vcf.gz"
        fix_merge_out="${sample_id}_merged_fixed_cipos_noTRA.vcf"
        '''
        echo !{vcf_DEL_cons} !{vcf_DUP_cons} !{vcf_TRA_cons} !{vcf_INV_cons}

        for type in DEL DUP TRA INV; do
                echo -e "\n ########## SORTING file consensus + FIX GENOTYPE + FIX CIPOS for: ${type} ########## \n"
		wc -l !{sample_id}_${type}_consensus.vcf
                echo -e "################################# ${type} \n --------SORT CONSENSUS--------"
                vcf-sort -c !{sample_id}_${type}_consensus.vcf > !{sample_id}_${type}_consensus_sorted.vcf
                wc -l !{sample_id}_${type}_consensus_sorted.vcf
		echo -e "################################# \n --------BGZIP--------"
                bgzip -f -c !{sample_id}_${type}_consensus_sorted.vcf > !{sample_id}_${type}_consensus_sorted.vcf.gz
                echo -e "################################# \n --------INDEX BCFTOOLS--------"
                bcftools index -t !{sample_id}_${type}_consensus_sorted.vcf.gz
                echo -e "################################# \n --------SURVIVOR > fixed--------"
                Rscript /opt/fix_SURVIVORgenotypes.R !{sample_id}_${type}_consensus_sorted.vcf.gz !{sample_id} !{sample_id}_${type}_consensus_sorted_fixed.vcf.gz
                echo -e "################################# \n --------GUNZIP--------"
                gunzip -f !{sample_id}_${type}_consensus_sorted_fixed.vcf.gz
                echo -e "################################# \n --------BGZIP--------"
                bgzip -f -c !{sample_id}_${type}_consensus_sorted_fixed.vcf > !{sample_id}_${type}_consensus_sorted_fixed.vcf.gz
                echo -e "################################# \n --------INDEX BCFTOOLS--------"
                bcftools index -t !{sample_id}_${type}_consensus_sorted_fixed.vcf.gz
                echo -e "################################# \n --------fix CIPOS--------"
                Rscript /opt/fix_CIPOS2.R !{sample_id}_${type}_consensus_sorted_fixed.vcf.gz !{sample_id}_${type}_consensus_sorted_fixed_cipos.vcf.gz
                gunzip -f !{sample_id}_${type}_consensus_sorted_fixed_cipos.vcf.gz
                bgzip -f -c !{sample_id}_${type}_consensus_sorted_fixed_cipos.vcf > !{sample_id}_${type}_consensus_sorted_fixed_cipos.vcf.gz
                bcftools index -t !{sample_id}_${type}_consensus_sorted_fixed_cipos.vcf.gz
        done

        echo -e "\n ########## BCFTOOLS CONCAT ON VCF with DEL DUP INV (no TRA because smoove can't genotype TRA)##########"
        #####Merge the 3 VCFs (no TRA because smoove genotype doesn't work for TRA) in one in order to later obtain one VCF for sample
        bcftools concat -a -O v -o !{fix_merge_out} \
                !{sample_id}_DEL_consensus_sorted_fixed_cipos.vcf.gz \
                !{sample_id}_DUP_consensus_sorted_fixed_cipos.vcf.gz \
                !{sample_id}_INV_consensus_sorted_fixed_cipos.vcf.gz
        '''
}
