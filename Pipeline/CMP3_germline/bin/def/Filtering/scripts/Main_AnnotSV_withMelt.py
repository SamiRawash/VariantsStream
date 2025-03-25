#! /usr/bin/env python3

import time, argparse, pandas, xlsxwriter, os
from manta_germ.fabioScript.provaFilteringannotsvvcf_withMelt import AnnotSV_AnnotationFiltering,STANDARD_OUT_COLUMNS
from vcf_tools.BamManipulation import find_binsize
from manta_germ.manta_melt.analysis_manta_melt import get3tables

#from manta_germ.fabioScript.analysis_manta_melt import get3tables

'''
Created on Mar 08, 2022

@author: afant
'''

def getListGenes(BED):
    selectedGenes=[]
    with open(BED) as bedfile:
        for line in bedfile:
            if not line.startswith("#"):
                splitted=line.split("\t")
                if len(splitted) > 3:
                    selectedGenes.append(splitted[3].strip())
    return selectedGenes


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This program does filters on SV using TSV annotated of Manta and TSV annotated of Melt (optional). The TSV must have been annotated with AnnotSV v3.2.3 and must have been processed with SVAFotate. \
    The output is a unique excel file with 4 sheets: 1. longCNV, 2. shortCNV, 3. INV, 4. TRA. In addition there are sheets 5. INS only in Manta and 5. MEI if MELT file is provided; \
    only sheet 4. INS if MELT file is NOT provided. An additional filter can be applied providing a bed file with genes of interest. Plots with samplot are done automatically for DEL and DUP when the bed file with genes of interest is provided. Automatically adds values of pHaplo, pTriplo. Output folder is created if doesn't exist.")
    parser.add_argument('-col','--column',action='store',type=str,help="List of columns to be used.", required=False, default=None)

    ### INPUT FILE
    parser.add_argument('-i','--input',action='store',type=str,help="Path of the input file (must be the output of AnnotSV).", required=True, default=None)
    parser.add_argument('-m','--melt',action='store',type=str,help="TSV input of melt.", required=False, default=None)
    parser.add_argument('-b','--bam',action='store',type=str,help="Path of the bam file (must correspond to the output of AnnotSV).", required=False, default=None)
    parser.add_argument('-g','--bed',action='store',type=str,help="Path of the bed file with list of coordinates of the genes of interest.", required=False, default=None)
    parser.add_argument('-d','--dosageFile',action='store',type=str,help="Path of the file with dosage scores (pHaplo, pTriplo)", required=False, default="/hpcshare/genomics/bioinfo_parabricks/BEDfiles/GeneDosage_converted_V43.tsv" )
    parser.add_argument('-c','--coverage',action='store',type=int,help="Coverage. Default is 30x ", required=False, default=30)
    parser.add_argument('-f','--frequency',action='store_true',help="Use this option to filter only for SV with AF < 0.01 ", required=False)

    ### OUTPUT - File
    parser.add_argument('-o','--outputfile',action='store',type=str,help="Path of the output folder where to put the output files.", required=True, default=None)
    arg=parser.parse_args()

    INPUT=arg.input
    #BAM=None
    BAM=arg.bam
    MELT=arg.melt
    FILE_DOSAGE=arg.dosageFile
    FREQ=arg.frequency
    OUTPUT_FOLDER=str(arg.outputfile)
    isExists = os.path.exists(OUTPUT_FOLDER)
    if not isExists:
        os.makedirs(OUTPUT_FOLDER)

    ### IF BED file is IDENTIFIED Load LIST GENE - ELSE FILTER on GENE NAMES is NOT APPLIED
    selectedGenes=getListGenes(arg.bed) if arg.bed is not None else None

    ### SELECT list columns
    COLUMNS=arg.column.split(',') if arg.column is not None else STANDARD_OUT_COLUMNS

    ### CREATE READER and LOAD DATA
    reader=AnnotSV_AnnotationFiltering(columns=COLUMNS)
    reader.verbose=True
    reader.target_coverage = arg.coverage
    #reader.meltfile = arg.melt
    reader.load(INPUT)
    reader.loadmanta(INPUT)
    reader.loadmelt(MELT)

    print("NEW ANALYSIS !!!!!!!!!!!!!!!!!!!!!!!!!!! \n \n")
    print(f"SAMPLE NAME: {reader.sample_name}")

    #nomi delle colonne di interesse, campi dell'INFO di cui voglio il valore, campi del FORMAT di cui voglio il valore
    #reader.unpackInfoFormat(reader.noImprecise[reader.dumped_columns], ["END"], ["GT", "PR", "SR"] )

    ### CREATE OUTPUT NAMES
    tsv_splitted = INPUT.split("/")
    filename = tsv_splitted[-1]
    #samplename = filename.replace("mantaGerm_annotsv.tsv", "")
    samplename = reader.sample_name

    #############################################################START - EXCEL OUTPUT#############################################################################

    OUT_XLSX = OUTPUT_FOLDER + "/" + samplename + "_StructuralVariants.xlsx"
    writer = pandas.ExcelWriter(OUT_XLSX, engine = 'xlsxwriter')
    print("OUTPUT_FOLDER: ", OUTPUT_FOLDER)
    print("EXCEL: ", OUT_XLSX)

    #SHORT
    print("##############################################################################\nSHORT CNV\n##############################################################################")
    sheet_short = "shortCNV"
    shortCNV, numberSVshort = reader.selectShortCNV(selectedGenes,folder_output=OUTPUT_FOLDER, bam_file=BAM, dosage=FILE_DOSAGE, filterAF=FREQ)
    shortCNV.to_excel(writer, sheet_name=sheet_short ,index =False)

    #LONG
    print("##############################################################################\nLONG CNV\n##############################################################################")
    sheet_long =  "longCNV"
    longCNV, numberSVlong = reader.selectLongCNV(selectedGenes, folder_output=OUTPUT_FOLDER, bam_file=BAM, dosage=FILE_DOSAGE, filterAF=FREQ)
    longCNV.to_excel(writer, sheet_name=sheet_long ,index =False)


    if arg.melt is not None:
        #INS only Manta
        print("##############################################################################\nINS\n##############################################################################")
        sheet_ins = "INS only in Manta"
        INS, numberINS = reader.selectINSonlyManta(selectedGenes=selectedGenes,folder_output=OUTPUT_FOLDER, filterAF=FREQ)
        INS.to_excel(writer, sheet_name=sheet_ins ,index =False)
        #MEI + INS
        print("##############################################################################\nMEI\n##############################################################################")
        sheet_mei = "MEI"
        MEI, numbersmei = reader.selectMelt(selectedGenes=selectedGenes,folder_output=OUTPUT_FOLDER, filterAF=FREQ)
        MEI.to_excel(writer, sheet_name=sheet_mei ,index =False)
    else:
        sheet_insNoMelt = "INS"
        INSwhenNoMelt, numberINS_whenNoMelt = reader.selectINSifNOMelt(selectedGenes=selectedGenes,folder_output=OUTPUT_FOLDER, filterAF=FREQ)
        INSwhenNoMelt.to_excel(writer, sheet_name=sheet_insNoMelt ,index =False)


    #INV
    print("##############################################################################\nINV\n##############################################################################")
    sheet_inv = "INV"
    INV, numbersinv = reader.selectINV(selectedGenes=selectedGenes,folder_output=OUTPUT_FOLDER, filterAF=FREQ)
    INV.to_excel(writer, sheet_name=sheet_inv ,index =False)

    #TRA
    print("##############################################################################\nTRA\n##############################################################################")
    sheet_tra = "TRA"
    TRA, numberstra = reader.selectTRA(selectedGenes=selectedGenes,folder_output=OUTPUT_FOLDER, filterAF=FREQ)
    TRA.to_excel(writer, sheet_name=sheet_tra ,index =False)



    ######
    print("GENERAL RESULTS: \n")
    print(f"Result {sheet_short}: ", numberSVshort)
    print(f"Result {sheet_long}: ", numberSVlong)
    print(f"Result {sheet_inv} : ", numbersinv)
    print(f"Result {sheet_tra} : ", numberstra)
    if arg.melt is not None:
        print(f"Result {sheet_ins}: ", numberINS)
        print(f"Result {sheet_mei} : ", numbersmei)
    else:
        print(f"MELT file NOT provided in the analysis \nResult {sheet_insNoMelt} : ", numberINS_whenNoMelt)


    writer.close()
    #############################################################END - EXCEL OUTPUT#############################################################################
    """
    #############################################################START - CSV OUTPUT#############################################################################
    print("OUTPUT_FOLDER: ", OUTPUT_FOLDER)
    #SHORT
    print("##############################################################################\nSHORT SV\n##############################################################################")
    output_name_shortCNV = OUTPUT_FOLDER + samplename + "manta_short.tsv"
    shortCNV, numberSVshort = reader.selectShortVariants(selectedGenes,folder_output=OUTPUT_FOLDER, bam_file=BAM)
    shortCNV.to_csv(output_name_shortCNV, sep = "\t" ,index =False)

    #LONG
    print("##############################################################################\nLONG SV\n##############################################################################")
    output_name_longCNV = OUTPUT_FOLDER + samplename + "manta_long.tsv"
    longCNV, numberSVlong = reader.selectLongVariants(selectedGenes, folder_output=OUTPUT_FOLDER, bam_file=BAM)
    #longCNV.to_csv(output_name_longCNV, sep = "\t", index = False)
    longCNV.to_csv(output_name_longCNV, sep = "\t" ,index =False)
    print("GENERAL RESULTS: \n")
    print(f"Result [Short] {output_name_shortCNV}: ", numberSVshort)
    print(f"Result [Long] {output_name_longCNV}: ", numberSVlong)
    """
    #############################################################END - CSV OUTPUT#############################################################################
