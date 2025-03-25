#! /usr/bin/env python3
import pandas, argparse, xlsxwriter, os
from manta_germ.manta_melt.analysis_manta_melt import get3tables, writeExcel
#from manta_germ.fabioScript.provaFilteringannotsvvcf_withMelt import takeID
#from manta_germ.fabioScript.Main_AnnotSV import getListGenes
list_mei=['ALU', 'HERVK', 'LINE1', 'SVA']
list_interesting = ['ALU', 'HERVK', 'LINE1', 'SVA', 'INS']

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
    parser = argparse.ArgumentParser(description="This script takes in input a tsv of Manta and a tsv of Melt and, optionally, a bed file with genes of interest. Produces a _Melt_Manta.xlsx with 3 sheets: one with MEI found both in Manta and Melt, one with INS only in Manta and one with MEI only in Melt.")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="TSV input of manta.", required=True, default=None)
    parser.add_argument('-m','--melt',action='store',type=str,help="TSV input of melt.", required=True, default=None)
    parser.add_argument('-o','--output',action='store',type=str,help="Output folder where I want to put my files.", required=True, default=None)
    parser.add_argument('-b','--bed',action='store',type=str,help="Bed files with the genes of interest.", required=False, default=None)
    arg=parser.parse_args()


    MELT=arg.melt
    MANTA=arg.input
    folder_output=str(arg.output)
    isExists = os.path.exists(folder_output)
    if not isExists:
        os.makedirs(folder_output)
    selectedGenes=getListGenes(arg.bed) if arg.bed is not None else None
    #print(selectedGenes)
    tsv_splitted = MANTA.split("/")
    filename = tsv_splitted[-1]
    samplename = filename.replace("mantaGerm_annotsv.tsv", "")
    out_file = folder_output + "/" + samplename + "Melt_Manta.xlsx"
    print("SAMPLE: ", samplename)

    #onlyIN_Manta, onlyIN_Melt, df = get3tables(tsv_Manta = MANTA, tsv_Melt=MELT, GENES=GENES)
    both, onlyIN_Manta, onlyIN_Melt, tabel_melt = get3tables(GENES=selectedGenes, tsv_Manta = MANTA, tsv_Melt=MELT)
    writeExcel(Manta=onlyIN_Manta, Mei=onlyIN_Melt, Both=both, out_file=out_file)
