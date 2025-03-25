process Annovar_Annotation {
input:
        tuple val(sample_id), path(anno_in)
        val reference_genome
        val bed_file
        val omim_file
output:
        tuple val(sample_id),path("${output}.vcf")
        tuple val(sample_id),path("${output}.txt")
script:
        sample_ID=anno_in.getSimpleName()
        THREADS=params.ANNOVAR_CPUS
        output="${sample_ID}_annovar.${reference_genome}_multianno"
        """
        perl /usr/local/bin/annovar/table_annovar.pl \
		"${anno_in}" \
		/home/dbs/ \
		-buildver "${reference_genome}" \
		-out "${sample_ID}_annovar" \
		-remove \
		-protocol refGene_latest,cytoBand,avsnp_latest,dbnsfp_latest,gnomad_genome_latest,icgc_latest,clinvar_latest,intervar_latest,dbscsnv_latest,dbnsfp_interpro_latest,bed,bed,bed \
		-operation g,r,f,f,f,f,f,f,f,f,r,r,r \
		-bedfile "${bed_file}","${omim_file}","${omim_file}" \
		-nastring . \
		-arg '-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,-hgvs,\
		--colsWanted 4, --colsWanted 13, --colsWanted 6' \
		-vcfinput \
		-polish \
		--thread ${THREADS}
        """
}
