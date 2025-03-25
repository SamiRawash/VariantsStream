#! /usr/bin/env python3
import pandas, os, numpy
from vcf_tools.BamManipulation import BAM_Manager, find_binsize
from manta_germ.manta_melt.analysis_manta_melt import get3tables


STANDARD_OUT_COLUMNS_oldAnnotSV=['AnnotSV_ID', 'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type', 'Samples_ID', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Replace_with_sample_name', 'Annotation_mode', 'Gene_name', 'Gene_count', 'Tx', 'Overlapped_CDS_percent','Exon_count','Location', 'Location2', 'Dist_nearest_SS', 'Nearest_SS_type','P_gain_source','P_loss_source','P_ins_source','P_snvindel_nb','P_snvindel_phen','B_gain_source','B_loss_source','B_ins_source','B_inv_source','TAD_coordinate', 'ACMG_class',  'ACMG', 'HI', 'TS',  'DDD_HI_percent', 'DDD_status', 'DDD_mode', 'DDD_consequence', 'OMIM_ID', 'OMIM_phenotype','OMIM_inheritance', 'OMIM_morbid', 'OMIM_morbid_candidate','GnomAD_pLI', 'ExAC_pLI', 'AnnotSV_ranking_score','AnnotSV_ranking_criteria']
STANDARD_OUT_COLUMNS=['AnnotSV_ID', 'SV_chrom', 'SV_start', 'SV_end', 'SV_length', 'SV_type', 'Samples_ID', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'Replace_with_sample_name', 'Annotation_mode', 'Gene_name', 'Gene_count', 'Tx', 'Overlapped_CDS_percent','Exon_count','Location', 'Location2', 'Dist_nearest_SS', 'Nearest_SS_type','P_gain_source','P_gain_phen','P_loss_source','P_loss_phen','P_ins_source','P_ins_phen','P_snvindel_nb','P_snvindel_phen','B_gain_source','B_gain_AFmax','B_loss_source','B_loss_AFmax','B_ins_source','B_ins_AFmax','B_inv_source','B_inv_AFmax','TAD_coordinate', 'ACMG_class',  'ACMG', 'HI', 'TS',  'DDD_HI_percent', 'GenCC_disease','GenCC_moi','GenCC_pmid','OMIM_ID', 'OMIM_phenotype','OMIM_inheritance', 'OMIM_morbid', 'OMIM_morbid_candidate','GnomAD_pLI', 'ExAC_pLI', 'AnnotSV_ranking_score','AnnotSV_ranking_criteria']

#bed_exon="/home/afant/Desktop/Parkinson/CLUSTERVDA/my_database/hg38_UCSC_exons_mergedName.bed.gz"
bed_exon="/home/afant/Desktop/Parkinson/CLUSTERVDA/my_database/202302_UCSC_CodingExonsPM100Hg38.bed.gz"
bed_all="/home/afant/Desktop/Parkinson/CLUSTERVDA/my_database/hg38_total_genesName.bed.gz"

class AnnotSV_AnnotationFiltering():

    def __init__(self, quality_pasing="PASS", imprecise="IMPRECISE", columns=STANDARD_OUT_COLUMNS):
        super().__init__()
        self.quality_passing=quality_pasing
        self.removed_rows=imprecise
        self.dumped_columns=columns
        #self.score_clingen=score_clingen
        self.longSV_threshold=500000
        #self.path_container_samplot="/home/afant/Desktop/Parkinson/container/samplot_latest.sif"
        self.path_container_samplot="samplot"

        ### DATA
        self.sample_name=None
        self.noImprecise=None
        self.verbose = False
        self.target_coverage = 30
        self.meltfile = None
        self.mantafile = None


    def __str__(self):
        return "CIAO"

    def splitTupla(self, tupla):
        id_sv = tupla[1]
        chrom=str(tupla[2])
        start=tupla[3]
        end=tupla[4]
        length=tupla[5]
        typ=tupla[6]
        chrom=chrom if len(chrom)>3 and chrom[0:3]=="chr" else f"chr{chrom}"
        return id_sv, chrom, start, end, length, typ

    def samplot_dump_sv(self, title, bam_file, chrom, start, end, typ, folder_output="", check_name=True, bed=None):
        bam_filename=bam_file.split("/")[-1]
        if check_name and not ( self.sample_name in bam_filename ) :
            raise Exception("Check of the name has fail!!!")

        out_png_samplot=folder_output + self.sample_name+"_"+title+".png"
        #cmd_str=f"singularity run {self.path_container_samplot} plot -n {self.sample_name} -c {chrom} -s {start} -e {end} -t {typ} -o {out_png_samplot} -b {bam_file} -A {bed}"
        #+ (f" -A {rmsk}" if int(end)-int(start)<2000000 else "" )
        cmd_str=f"{self.path_container_samplot} plot -n {self.sample_name} -c {chrom} -s {start} -e {end} -t {typ} -o {out_png_samplot} -b {bam_file} -A {bed}"
        print(cmd_str)
        os.system(cmd_str)

    def load(self, file, separator='\t'):
        f = pandas.read_csv(file, sep=separator)
        self.sample_name=f.columns[14]

        ### GET COLUMN FORMAT names
        self.dumped_columns=[ self.sample_name if x=='Replace_with_sample_name' else x for x in STANDARD_OUT_COLUMNS ]

        f=f[self.dumped_columns]
        self.noImprecise_pre = f.loc[(~f['INFO'].str.contains(self.removed_rows) & (f['FILTER'] == self.quality_passing))]
        self.noImprecise = self.unpackInfoFormat(self.noImprecise_pre, ["END","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])

    def loadmelt(self, file):
        self.meltfile = file
    def loadmanta(self, file):
        self.mantafile = file

    def unpackInfoFormat(self, table, columnListINFO=[], columnListFORMAT=[]):
        new_columns={}
        for id in columnListINFO:
            new_columns[id]=[]
        for id in columnListFORMAT:
            new_columns[id]=[]
        for row in table[['INFO','FORMAT',self.sample_name]].itertuples():
            info=self.parseInfo(row[1])
            format=self.parseFormat(row[2],row[3])
            for id in columnListINFO :
                new_columns[id].append( info[id] if id in info else None )
            for id in columnListFORMAT :
                new_columns[id].append( format[id] if id in format else None )
        counter=15
        df=pandas.DataFrame(table)
        for id in columnListINFO:
            df.insert(counter,id,new_columns[id],True)
            counter+=1
        for id in columnListFORMAT:
            df.insert(counter,id,new_columns[id],True)
            counter+=1
        return df

    def chose_bin(self,x):
        if x<=100 :
            return 3
        elif x<=1000 :
            return 2
        return 1


    def produceHist_Table(self, bam_file=None, table=None, folder_output="", bedfile=None, selectedGenes=None, dosage=None):
        dictionary_validation_hist={}
        if bam_file is not None :
            BAM=BAM_Manager(bam_file)
            BAM.verbose = self.verbose
            for tupla in table[["AnnotSV_ID","SV_chrom","SV_start","SV_end","SV_length","SV_type"]].itertuples() :
                id_sv, chrom, start, end, length, typ = self.splitTupla(tupla)
                if (typ == "DEL" or typ == "DUP") and type(selectedGenes) == list:
                    print("GENERATE SAMPLOT PLOT: ", id_sv)
                    self.samplot_dump_sv(id_sv, bam_file, chrom, start, int(start+abs(length)), typ, folder_output=folder_output, bed=bedfile)
        print("col no Imprecise: ", self.noImprecise.columns)
        tabellaconID_pre = self.takeID(None, None, table, self.noImprecise)
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!tabellaconID_pre \n ", tabellaconID_pre.columns)
        tabellaconID=pandas.DataFrame(tabellaconID_pre, index=None)

        tabellaconID["GENE_INTEREST"] = ""
        #Se ho bed con geni
        if selectedGenes is not None:
            for gene in tabellaconID['Gene_name']:
                tabellaconID.loc[((tabellaconID['Gene_name'] == gene) & (tabellaconID['Annotation_mode'] == "split") & (selectedGenes.count(gene) >0)), 'GENE_INTEREST'] = "*"
        if dosage is not None:
            table_geneDosage = pandas.read_csv(dosage, sep="\t")
            tabellaconID_new = pandas.merge(tabellaconID, table_geneDosage, on=['Gene_name'], how='left')
        return tabellaconID_new


    def selectShortCNV(self, selectedGenes, folder_output="", bam_file=None, dosage=None, filterAF=False) :
        #saved_args = locals()
        if selectedGenes is not None :
            #FILTER ON SV_LENGTH < 500.000 and SV in GENE
            shortCNV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & ((self.noImprecise['SV_type'] == "DEL") | (self.noImprecise['SV_type'] == "DUP")) & (self.noImprecise['SV_length'].abs().le(self.longSV_threshold) ) & ([ len( set(str(x).split(";")).intersection(selectedGenes))>0 for x in self.noImprecise['Gene_name'] ] )]
        else :
            #FILTER ON SV_LENGTH < 500.000
            shortCNV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full")& ((self.noImprecise['SV_type'] == "DEL") | (self.noImprecise['SV_type'] == "DUP")) &  (self.noImprecise['SV_length'].abs().le(self.longSV_threshold) )]

        #Seleziono la tabella con parsate le colonne di interesse da INFO e FORMAT
        #table_parsed=self.unpackInfoFormat(shortCNV, ["END","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])
        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = shortCNV
        elif filterAF==True:
            table_parsed_filteredAF = shortCNV.loc[shortCNV['Max_AF'].astype(float) < 0.01]

        tabellaconID_short = self.produceHist_Table( bam_file=bam_file, table=table_parsed_filteredAF, folder_output=folder_output, bedfile=bed_exon, selectedGenes=selectedGenes, dosage=dosage)
        #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n !!!!!!!!!!!!!!!!!!!!!!tabellaconID_short: ", tabellaconID_short.columns)

        #print("tabellaconID_short", tabellaconID_short)
        #Number of SVs in genelist
        numbersSV = tabellaconID_short['AnnotSV_ID'].nunique()
        return tabellaconID_short, numbersSV


    def selectLongCNV(self, selectedGenes, folder_output="", bam_file=None, dosage=None, filterAF=False):
        if selectedGenes is not None:
            #FILTER ON SV_LENGTH > 500.000
            longCNV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & ((self.noImprecise['SV_type'] == "DEL") | (self.noImprecise['SV_type'] == "DUP")) & (self.noImprecise['SV_length'].abs().gt(self.longSV_threshold)) & ([ len( set(str(x).split(";")).intersection(selectedGenes))>0 for x in self.noImprecise['Gene_name'] ] ) ]
        else:
            longCNV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & ((self.noImprecise['SV_type'] == "DEL") | (self.noImprecise['SV_type'] == "DUP")) & (self.noImprecise['SV_length'].abs().gt(self.longSV_threshold) )]

        #table_parsed=self.unpackInfoFormat(longCNV, ["END","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])
        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = longCNV
        elif filterAF==True:
            table_parsed_filteredAF = longCNV.loc[longCNV['Max_AF'].astype(float) < 0.01]

        tabellaconID_long = self.produceHist_Table( bam_file=bam_file, table=table_parsed_filteredAF, folder_output=folder_output, bedfile=bed_all, selectedGenes=selectedGenes, dosage=dosage)

        numbersSV = tabellaconID_long['AnnotSV_ID'].nunique()
        return tabellaconID_long, numbersSV


    def selectINV(self, selectedGenes, folder_output="", filterAF=False):
        if selectedGenes is not None:
            #FILTER ON SV_LENGTH > 500.000
            INV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "INV")  & ([ len( set(str(x).split(";")).intersection(selectedGenes))>0 for x in self.noImprecise['Gene_name'] ] ) ]
        else:
            INV = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "INV")]
        tabellaconID_pre = self.takeID(None, None, INV, self.noImprecise)
        tabellaconID=pandas.DataFrame(tabellaconID_pre, index=None)
        tabellaconID["GENE_INTEREST"] = ""
        if selectedGenes is not None:
            for gene in tabellaconID['Gene_name']:
                tabellaconID.loc[((tabellaconID['Gene_name'] == gene) & (tabellaconID['Annotation_mode'] == "split") & (selectedGenes.count(gene) >0)), 'GENE_INTEREST'] = "*"

        #table_parsed = self.unpackInfoFormat(tabellaconID, ["END","CHR2","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])

        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = tabellaconID
        elif filterAF==True:
            table_parsed_filteredAF = tabellaconID.loc[tabellaconID['Max_AF'].astype(float) < 0.01]
        numbersSV = table_parsed_filteredAF['AnnotSV_ID'].nunique()
        return table_parsed_filteredAF, numbersSV

    def selectTRA(self, selectedGenes, folder_output="", filterAF=False):
        if selectedGenes is not None:
            #print(("N[chr" + str(self.noImprecise['SV_chrom']) + ":1[" ))
            TRA = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "TRA") & (~self.noImprecise['ALT'].str.endswith(":1[")) & ([ len( set(str(x).split(";")).intersection(selectedGenes))>0 for x in self.noImprecise['Gene_name'] ] ) ]
        else:
            #print(self.noImprecise.loc[("N[chr" + str(self.noImprecise['SV_chrom']) + ":1[")])
            TRA = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "TRA") & (~self.noImprecise['ALT'].str.endswith(":1["))]
        tabellaconID_pre = self.takeID(None, None, TRA, self.noImprecise)
        tabellaconID=pandas.DataFrame(tabellaconID_pre, index=None)
        tabellaconID["GENE_INTEREST"] = ""
        if selectedGenes is not None:
            for gene in tabellaconID['Gene_name']:
                tabellaconID.loc[((tabellaconID['Gene_name'] == gene) & (tabellaconID['Annotation_mode'] == "split") & (selectedGenes.count(gene) >0)), 'GENE_INTEREST'] = "*"

        #table_parsed = self.unpackInfoFormat(tabellaconID, ["END","CHR2","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])

        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = tabellaconID
        elif filterAF==True:
            table_parsed_filteredAF = tabellaconID.loc[tabellaconID['Max_AF'].astype(float) < 0.01]
        numbersSV = table_parsed_filteredAF['AnnotSV_ID'].nunique()
        return table_parsed_filteredAF, numbersSV

    def selectINSifNOMelt(self, selectedGenes=None, folder_output="", filterAF=False):
        if selectedGenes is not None:
            #FILTER ON SV_LENGTH > 500.000
            INS_noMelt = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "INS") & ([ len( set(str(x).split(";")).intersection(selectedGenes))>0 for x in self.noImprecise['Gene_name'] ] ) ]
        else:
            INS_noMelt = self.noImprecise.loc[(self.noImprecise['Annotation_mode'] == "full") & (self.noImprecise['SV_type'] == "INS")]
        tabellaconID_pre = self.takeID(None, None, INS_noMelt, self.noImprecise)
        tabellaconID=pandas.DataFrame(tabellaconID_pre, index=None)
        tabellaconID["GENE_INTEREST"] = ""
        if selectedGenes is not None:
            for gene in tabellaconID['Gene_name']:
                tabellaconID.loc[((tabellaconID['Gene_name'] == gene) & (tabellaconID['Annotation_mode'] == "split") & (selectedGenes.count(gene) >0)), 'GENE_INTEREST'] = "*"

        #table_parsed = self.unpackInfoFormat(tabellaconID, ["END","CHR2","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])
        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = tabellaconID
        elif filterAF==True:
            table_parsed_filteredAF = tabellaconID.loc[tabellaconID['Max_AF'].astype(float) < 0.01]

        numbersSV = table_parsed_filteredAF['AnnotSV_ID'].nunique()
        return table_parsed_filteredAF, numbersSV

    def selectINSonlyManta(self, selectedGenes=None, folder_output="", filterAF=False):
        both, onlyIN_Manta_pre, onlyIN_Melt, table_melt  = get3tables(tsv_Manta=self.mantafile, tsv_Melt=self.meltfile, GENES=selectedGenes)
        #table_parsed = self.unpackInfoFormat(onlyIN_Manta, ["END","CHR2","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])
        onlyIN_Manta = self.unpackInfoFormat(onlyIN_Manta_pre, ["END","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])
        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = onlyIN_Manta
        elif filterAF==True:
            table_parsed_filteredAF = onlyIN_Manta.loc[onlyIN_Manta['Max_AF'].astype(float) < 0.01]
        numbersSV = table_parsed_filteredAF['AnnotSV_ID'].nunique()
        return table_parsed_filteredAF, numbersSV


    def selectMelt(self, selectedGenes=None, folder_output="", filterAF=False):
        both, onlyIN_Manta, onlyIN_Melt, table_melt_pre = get3tables(tsv_Manta=self.mantafile, tsv_Melt=self.meltfile, GENES=selectedGenes)
        table_melt = self.unpackInfoFormat(table_melt_pre, ["END","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP"],  ["GT", "DP", "RP", "AP", "RS", "AS", "AB", "DHGT" ])

        tabellamelt_complete_pre=pandas.DataFrame(table_melt, index=None)
        tabellamelt_complete_pre["GENE_INTEREST"] = ""
        tabellamelt_complete_pre["Manta_confirmation"] = ""
        if filterAF==False:
            tabellamelt_complete = tabellamelt_complete_pre
        elif filterAF==True:
            tabellamelt_complete = tabellamelt_complete_pre.loc[tabellamelt_complete_pre['Max_AF'].astype(float) < 0.01]

        ID_both = both[both.columns[0]]
        LIST_ID_both = list(set(ID_both))
        for myid in tabellamelt_complete['AnnotSV_ID']:
            tabellamelt_complete.loc[tabellamelt_complete['AnnotSV_ID'].isin(LIST_ID_both), 'Manta_confirmation'] = "yes"
        if selectedGenes is not None:
            for gene in tabellamelt_complete['Gene_name']:
                tabellamelt_complete.loc[(tabellamelt_complete['Gene_name'] == gene) & (tabellamelt_complete['Annotation_mode'] == "split") & (selectedGenes.count(gene) >0), 'GENE_INTEREST' ] = "*"

        #table_parsed = self.unpackInfoFormat(tabellamelt_complete, ["END","CHR2","Max_AF", "Max_Het", "Max_HomAlt", "Best_1000G_ID", "Best_1000G_OFP", "Best_gnomAD_ID", "Best_gnomAD_OFP", "Best_CCDG_ID", "Best_CCDG_OFP", "ASSESS", "RP", "SR"],  ["GT", "DP", "AD"])
        #Se viene data l'opzione -f filtra per AF < 0.01
        if filterAF==False:
            table_parsed_filteredAF = tabellamelt_complete
        elif filterAF==True:
            table_parsed_filteredAF = tabellamelt_complete.loc[tabellamelt_complete['Max_AF'].astype(float) < 0.01]

        numbersmei = table_parsed_filteredAF['AnnotSV_ID'].nunique()
        return table_parsed_filteredAF, numbersmei

    #By default takes a table (shortTable), an sv_type (e.g. DEL) and produced a table with filters in filter_sv () applied
    def takeID(self, sv_type, filter_sv, shortTable, completetable): #, full_table):
        column_id = None
        set_values = None
        if filter_sv is not None :
          column_id = filter_sv[0]
          set_values = filter_sv[1]
        ### Filter for the SV TYPE of INTEREST with a certain SCORE
        if sv_type is None:
          tmp_table = shortTable
        elif set_values is not None:
          #tmp_table=shortTable.loc[(shortTable['SV_type'] == sv_type) & (shortTable[column_id].isin(set_values))]
          SVtable = shortTable.loc[(shortTable['SV_type'] == sv_type)]
          IDofInterest = SVtable[SVtable.columns[0]]
          #print(DEL_ID)
          tableWithSelectedID = completetable.loc[completetable[completetable.columns[0]].isin(IDofInterest)]
          tmp_table = tableWithSelectedID.loc[(tableWithSelectedID['Annotation_mode'] == "split") & (tableWithSelectedID[column_id].isin(set_values)) ]
          #print(tmp_table)
          #DEL_OK = tableWithSelectedID_filtered[tableWithSelectedID_filtered.columns[0]]
          #tmp_table = self.noImprecise.loc[self.noImprecise[]]
        else :
          tmp_table = shortTable.loc[(shortTable['SV_type'] == sv_type)]

        ### COLLECT the ID of the SV of interest
        DEL_ID = tmp_table[tmp_table.columns[0]] #.itertuples()
        #return full_table.loc[full_table[full_table.columns[0]].isin(DEL_ID)]
        return completetable.loc[completetable[completetable.columns[0]].isin(DEL_ID)]


    #def takeID(self, SV_type_table, full_table):
    #    #DEL_ID=SV_type_table[colonne[0]]
    #    DEL_ID=SV_type_table[SV_type_table.columns[0]] #.itertuples()
    #    return full_table.loc[full_table[full_table.columns[0]].isin(DEL_ID)]

    def parseInfo(self, info_field):
        infoDict = {}
        each = info_field.split(';')
        for i in each:
            splitted= i.split("=")
            key= splitted[0]
            value=splitted[1] if len(splitted) > 1 else True
            infoDict[key] = value
        return infoDict


    def parseFormat(self, column_key,column_value):
        keys=column_key.strip().split(":")
        values=column_value.strip().split(":")
        formatDict = {}
        for n,key in enumerate(keys):
            formatDict[key.strip()]=values[n]
        return formatDict
