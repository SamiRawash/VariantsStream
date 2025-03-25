#! /usr/local/bin/python3

import argparse
import sys

def split_line_by_separator(vcf_line,separator="\t"):
    splitted_row=vcf_line.split(separator)
    return splitted_row

def get_genes_from_vcf_line(vcf_line,starting_string):
    genes_from_vcf_line=[]
    info=split_line_by_separator(vcf_line)[7]
    splitted_info=split_line_by_separator(info,";")
    for element in splitted_info:
        if element.startswith(starting_string):
            splitted_ann=split_line_by_separator(element,",")
            for el in splitted_ann:
                splitted_ann_pipe=split_line_by_separator(el,"|")[3]
                genes_from_vcf_line.append(splitted_ann_pipe)
    return genes_from_vcf_line

def get_genes_from_file(genes_file,separator=",",column=0):
    genes_from_file=[]
    line=genes_file.readline().rstrip()
    row=1 
    while line != "":
        splitted_line=split_line_by_separator(line,separator)
        column_number=len(splitted_line)
        if column_number > column:
            gene=splitted_line[column]
            #print(gene)
            genes_from_file.append(gene)
        else:
            return -1,row 

        line=genes_file.readline().rstrip()
        row=row+1
    return genes_from_file
   
      

if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Get gene from VCF")
    parser.add_argument('-v','--vcf',action='store',type=str,help='Path to the VCF file',required=True,default=None)
    parser.add_argument('-g','--genes',action='store',type=str,help='File containig genes to be selected from the VCF',required=True,default=None)
    parser.add_argument('-s','--separator',action='store',type=str,help='Separator of columns in the genes file',required=False,default=",")
    parser.add_argument('-c','--column',action='store',type=int,help='Column selected from the genes file',required=False,default=1)
    parser.add_argument('-e','--exclude',action='store_true',help='This option exclude from the vcf the genes included --genes option',required=False)
    parser.add_argument('-a','--annotator',action='store',type=str,help='Annotator chooser, s for snpeff and v for vep (s default value)',required=False,default="s")
    args=parser.parse_args()

    vcf_file_path=args.vcf
    genes_file_path=args.genes
    separator_genes_file=args.separator
    column_genes_file=args.column -1  
    exclude_genes=args.exclude
    annotator =args.annotator

    annotator_dictionary={"v":"CSQ=","s":"ANN="}
    try:
        starting_string=annotator_dictionary[annotator]
    except:
        sys.exit("The annotator option {} doesn't exist".format(annotator))

    
    try:
       vcf_file=open(vcf_file_path)
    except:
       sys.exit("The file {} could not be found".format(vcf_file_path)) #.format si riferisce alle {} devo dirgli cosa mettere nelle {}
    try:
        genes_file=open(genes_file_path)
    except:
        sys.exit("The file {} could not be found".format(genes_file_path)) 

    genes_from_file=get_genes_from_file(genes_file,separator_genes_file,column_genes_file)
    
    if genes_from_file[0] == -1:
        sys.exit("The number of columns is smaller than the number you have passed with -c")

    genes_from_file_set=set(genes_from_file)


    line=vcf_file.readline()
    n=1  
    
    while line!="":
        if line[0]=="#":
           
           print(line.rstrip())  
        else:
            genes_from_vcf_line=get_genes_from_vcf_line(line,starting_string)
            genes_from_vcf_set=set(genes_from_vcf_line)   
            intersection = genes_from_vcf_set & genes_from_file_set  
           
                                                                                                                                                             
            if len(intersection) > 0 and exclude_genes==False:
                print(line.rstrip())
              
            elif len(intersection) == 0 and exclude_genes==True:
                print(line.rstrip())
                                        
        line=vcf_file.readline()

    

    vcf_file.close()
    genes_file.close()

