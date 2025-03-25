#! /usr/bin/env python3
import pandas, argparse, xlsxwriter, os
#from manta_germ.fabioScript.provaFilteringannotsvvcf_withMelt import takeID
list_mei=['ALU', 'HERVK', 'LINE1', 'SVA']
list_interesting = ['ALU', 'HERVK', 'LINE1', 'SVA', 'INS']
STANDARD_OUT_COLUMNS=['AnnotSV_ID', 'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type', 'Samples_ID', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Replace_with_sample_name', 'Annotation_mode', 'Gene_name', 'Gene_count', 'Tx', 'Overlapped_CDS_percent','Exon_count','Location', 'Location2', 'Dist_nearest_SS', 'Nearest_SS_type','P_gain_source','P_gain_phen','P_loss_source','P_loss_phen','P_ins_source','P_ins_phen','P_snvindel_nb','P_snvindel_phen','B_gain_source','B_gain_AFmax','B_loss_source','B_loss_AFmax','B_ins_source','B_ins_AFmax','B_inv_source','B_inv_AFmax','TAD_coordinate', 'ACMG_class',  'ACMG', 'HI', 'TS',  'DDD_HI_percent', 'GenCC_disease','GenCC_moi','GenCC_pmid','OMIM_ID', 'OMIM_phenotype','OMIM_inheritance', 'OMIM_morbid', 'OMIM_morbid_candidate','GnomAD_pLI', 'ExAC_pLI', 'AnnotSV_ranking_score','AnnotSV_ranking_criteria']


def get3tables(GENES=None, tsv_Manta = None, tsv_Melt=None):
    #Get a table with only the genes of interest if provided
    #table_all = pandas.read_csv(tableMergedMantaMelt, sep = "\t")
    table_manta_pre = pandas.read_csv(tsv_Manta, sep = "\t")
    table_melt_pre = pandas.read_csv(tsv_Melt, sep = "\t")
    sample_name_manta=table_manta_pre.columns[14]
    sample_name_melt=table_melt_pre.columns[14]
    print("THE SAMPLE NAME:", sample_name_manta)
    dumped_columns_manta=[ sample_name_manta if x=='Replace_with_sample_name' else x for x in STANDARD_OUT_COLUMNS ]
    dumped_columns_melt=[ sample_name_melt if x=='Replace_with_sample_name' else x for x in STANDARD_OUT_COLUMNS ]

    table_manta_pre=table_manta_pre[dumped_columns_manta]
    table_melt_pre=table_melt_pre[dumped_columns_melt]

    if GENES is None :
        table_manta = table_manta_pre.loc[ (~table_manta_pre['INFO'].str.contains("IMPRECISE"))   & (table_manta_pre['FILTER'] == 'PASS') & (table_manta_pre['SV_type'] == 'INS' )]
        table_manta_onlyfull = table_manta.loc[table_manta['Annotation_mode'] == "full"]
        table_melt = table_melt_pre.loc[ table_melt_pre['FILTER'] == 'PASS']
        table_melt_onlyfull = table_melt.loc[table_melt['Annotation_mode'] == "full"]
    else: #AGGIUNGI QUI I ONLY FULL
        table_manta = table_manta_pre.loc[ (~table_manta_pre['INFO'].str.contains("IMPRECISE"))   & (table_manta_pre['FILTER'] == 'PASS') & (table_manta_pre['SV_type'] == 'INS' ) & ([ len( set(str(x).split(";")).intersection(GENES))>0 for x in table_manta_pre['Gene_name'] ])]
        table_manta_onlyfull = table_manta.loc[table_manta['Annotation_mode'] == "full"]
        table_melt = table_melt_pre.loc[ (table_melt_pre['FILTER'] == 'PASS') & ([ len( set(str(x).split(";")).intersection(GENES))>0 for x in table_melt_pre['Gene_name'] ])]
        table_melt_onlyfull = table_melt.loc[table_melt['Annotation_mode'] == "full"]

    #ONLY MEI
    #tableMEI = table.loc[table["SV_type"].isin(list_mei) ]
    #print("TOTAL MEI in MELT: ", tableMEI['AnnotSV_ID'].nunique())
    print("TOTAL MEI in MELT: ", table_melt.loc[table_melt['Annotation_mode'] == "full"].shape[0])

    #TABLE WITH INS OF MANTA
    #table_INS = table.loc[(table['SV_type'] == 'INS' ) & (table['Annotation_mode'] == "full") ]
    #print("TOTAL INS in MANTA: ", table_INS['AnnotSV_ID'].nunique())
    print("TOTAL INS in MANTA: ", table_manta.loc[table_manta['Annotation_mode'] == "full"].shape[0])

    #empty dataframe to fill with intersection of ins manta + mei melt
    both = pandas.DataFrame(columns=STANDARD_OUT_COLUMNS,  dtype=object)

    for tupla in table_melt_onlyfull[['SV_chrom', 'SV_start']].itertuples():
        chrom=tupla[1]
        start=int(tupla[2])
        infRange=int(start -20)
        supRange=int(start +20)
        #THIS TABLE CONTAINS UNION MELT MANTA
        #INSuMEI = table_manta.loc[ (table_manta['Annotation_mode'] == "full") & (table_manta['SV_chrom'] == chrom) & (table_manta['SV_start'] <= supRange) & (table_manta['SV_start'] >= infRange)]
        #MEIuINS = table_melt_onlyfull.loc[(table_melt_onlyfull['SV_chrom'] == chrom) & (table_melt_onlyfull['SV_start'] == start)]
        INSuMEI = table_manta.loc[ (table_manta['SV_chrom'] == chrom) & (table_manta['SV_start'] <= supRange) & (table_manta['SV_start'] >= infRange)]
        MEIuINS = table_melt.loc[(table_melt['SV_chrom'] == chrom) & (table_melt['SV_start'] == start)]
        if len(INSuMEI) >= 1:
            both = pandas.concat([both, INSuMEI], ignore_index=True, sort=False)
            both = pandas.concat([both, MEIuINS], ignore_index=True, sort=False)
    #print((both.loc[both['SV_type'] == "INS", 'AnnotSV_ID']).nunique() )
    print("INS FOUND IN BOTH: ", ((both.loc[both['SV_type'] == "INS", 'AnnotSV_ID']).nunique() ) )
    print("MEI FOUND IN BOTH: ", ((both.loc[both['SV_type'] != "INS", 'AnnotSV_ID']).nunique() ) )

    ID_BOTH = both[both.columns[0]]
    LIST_ID_BOTH = list(set(ID_BOTH))
    #print("BOTH", LIST_ID_BOTH)

    #INS ONLY IN MANTA -OK
    onlyIN_Manta = table_manta.loc[~table_manta['AnnotSV_ID'].isin(LIST_ID_BOTH)]
    #onlyIN_Manta.to_csv("onlyIN_Manta.csv", sep = "\t" ,index =False)
    print("INS ONLY in MANTA:", onlyIN_Manta['AnnotSV_ID'].nunique() )

    #MEI ONLY IN MELT
    onlyIN_Melt = table_melt.loc[~table_melt['AnnotSV_ID'].isin(LIST_ID_BOTH)]
    #print("ONLY MELT:",list(set(onlyIN_Melt[onlyIN_Melt.columns[0]])))
    print("MEI ONLY in MELT:", onlyIN_Melt['AnnotSV_ID'].nunique() )

    #Percentage of MEI found by Manta
    if (int((table_melt.loc[table_melt['Annotation_mode'] == "full"].shape[0])) > 0) :
        the_percentage = (( ((both.loc[both['SV_type'] != "INS", 'AnnotSV_ID']).nunique() )  )*100 /  (table_melt.loc[table_melt['Annotation_mode'] == "full"].shape[0]) )
        print( format(the_percentage, ".2f"), "% of MEI of MELT found also by MANTA \n")
    else:
        print("There are no MEI found by Melt -- no intersection with Manta")

    #return both, onlyIN_Manta, onlyIN_Melt, table_melt_onlyfull, table_melt, table_manta
    return both, onlyIN_Manta, onlyIN_Melt, table_melt

def writeExcel(Manta=None, Mei=None, Both=None, out_file=None):
    writer = pandas.ExcelWriter(out_file, engine = 'xlsxwriter')
    Both.to_excel(writer, sheet_name='MEI_U_INS' ,index =False)
    #print(ID_mantaMelt)
    #out_manta = folder_output + samplename + "INS_only.tsv"
    #onlyIN_Manta.to_csv(out_manta, sep = "\t", index =False)
    Manta.to_excel(writer, sheet_name='INS_only' ,index =False)
    #out_melt = folder_output + samplename + "MEI_only.tsv"
    #onlyIN_Melt.to_csv(out_melt, sep = "\t", index =False)
    Mei.to_excel(writer, sheet_name='MEI_only' ,index =False)
    writer.save()
