#! /usr/bin/env python3
'''
Created on 28 Feb 2023

@author: flanduzzi
'''
DEFAULT_ROW_SEPARATOR="\n"
DEFAULT_COLUMN_SEPARATOR="\t"
DEFAULT_COLUMN_INDEX_GENE=4
DEFAULT_COLUMN_INDEX_TRANSCRIPT=3
DEFAULT_FILE_CANONICAL='/opt/lib_5000GenomiVdA/vcf_tools/data/UCSCGenomeBrowser_CanonicalFrom_kgXref_20230228.csv' 
## 'data/UCSCGenomeBrowser_CanonicalFrom_kgXref_20230228.csv'
## "data/UCSCGenomeBrowser_EnsRefSeq_202302241905.csv"


def loadCSV(source_path, separator=DEFAULT_COLUMN_SEPARATOR, column=0, mode="r", header=True):
    data={}
    with open(source_path, mode=mode) as csv:
        if header :
            csv.readline()
        for line in csv.readlines() :
            if line[0]!="#" and line.strip()!="" :
                splitted=line.strip().split(separator)
                key=splitted[column].strip()
                data[key]=splitted
    return data


def loadCanonicalGeneFile(source_path, separator=DEFAULT_COLUMN_SEPARATOR, column_gene=0, column_transcript=1, mode="r", header=True, gene_subset=None):
    canonical_list={}
    with open(source_path, mode=mode) as csv:
        if header :
            csv.readline()
        counter=0
        for line in csv.readlines() :
            if line[0]!="#" and line.strip()!="" :
                counter+=1
                splitted=line.strip().split(separator)
                key=splitted[column_gene].strip()
                if gene_subset is None or key in gene_subset :
                    if key in canonical_list :
                        print("DOPPIONE: ")
                        print(f" - {key} : {canonical_list[key]}")
                        print(f" - {key} : {splitted[column_transcript]}")
                    else:
                        canonical_list[key]=splitted[column_transcript].strip()
    print(f"Canonical Genes [{counter}]: {len(canonical_list)}")
    return canonical_list            
    
def execute_loading():
    parser = argparse.ArgumentParser(prog="GenerateCanonicalGeneList", description="Demo to test the program functionalities.")
    parser.add_argument('-cs','--canonical_transcript_file',action='store_true', help=f"Path of the file containing the canonical transcript [DEFAULT: {DEFAULT_FILE_CANONICAL}].", default=DEFAULT_FILE_CANONICAL)
    parser.add_argument('-cg','--column_gene',action='store',type=int, help=f"Index of the column containing the GENE names convention is to start column count from zero [DEFAULT: {DEFAULT_COLUMN_INDEX_GENE}].", default=DEFAULT_COLUMN_INDEX_GENE)
    parser.add_argument('-ct','--column_transcript',action='store',type=int, help=f"Index of the column containing the TRANSCRIPT names convention is to start column count from zero [DEFAULT: {DEFAULT_COLUMN_INDEX_TRANSCRIPT}].", default=DEFAULT_COLUMN_INDEX_TRANSCRIPT)
   
    parser.add_argument('-i','--gene_list',action='store',type=str, help="Path of the file containing the canonical transcript [DEFAULT: None].", default=None)
    parser.add_argument('-s','--separator',action='store',type=str,help="Column separator in the gene_list file [DEFAULT: '{:}'].".format(DEFAULT_COLUMN_SEPARATOR), required=False, default=DEFAULT_COLUMN_SEPARATOR)
    parser.add_argument('-c','--column',action='store',type=int, help=f"Index of the column containing the GENE OF INTEREST names convention is to start column count from zero [DEFAULT: {DEFAULT_COLUMN_INDEX_GENE}].", default=DEFAULT_COLUMN_INDEX_GENE)
    parser.add_argument('-noheader','--noheader',action='store_true', help="GeneList file contains header [DEFAULT: False].", default=False)
    
    parser.add_argument('-o','--output',action='store',type=str, help="Path of the output file [DEFAULT: None].", default=None)
    arg=parser.parse_args()
    
    geneNameList=None
    if arg.gene_list is not None and arg.gene_list != "" :
        geneNameList=loadCSV(arg.gene_list, arg.separator, column=arg.column, header=(not arg.noheader))
    else :
        print("No specific gene list has been loaded all canonical genes will be used.")
    canonical_subset=loadCanonicalGeneFile(arg.canonical_transcript_file, column_gene=arg.column_gene, column_transcript=arg.column_transcript, gene_subset=None if geneNameList is None else set(geneNameList.keys()))
    
    print(f"Total Canonical transcript in the list: {len(canonical_subset)}")
    if geneNameList is not None :
        requested_gene=set(geneNameList.keys())
        find_genes=set(canonical_subset.keys())
        missed=list(requested_gene-find_genes)
        missed.sort()
        if len(missed)>0 :
            print(f"In the given list [{len(geneNameList)}] there are {len(missed)} not found: "+", ".join(missed))
    
    keys=list(canonical_subset.keys())
    keys.sort()
    if arg.output is not None :
        with open(arg.output, "w") as outfile :
            outfile.write(f"#Gene\tCanonicalTranscript\n")
            for el in keys :
                outfile.write(f"{el}\t{canonical_subset[el]}\n")            
    else :
        for el in keys :
            print(f"{el}\t{canonical_subset[el]}")
            
    return canonical_subset
                
if(__name__=='__main__'):
    print(("#"*54)+"\n### TEST GenerateCanonicalGeneList Functionalities ###\n"+("#"*54))
    import argparse
    geneNameList=execute_loading()
    
    from collections import deque
    shortList=deque([],maxlen=10)
    for el in geneNameList :
        shortList.append(el)
        #print(shortList)