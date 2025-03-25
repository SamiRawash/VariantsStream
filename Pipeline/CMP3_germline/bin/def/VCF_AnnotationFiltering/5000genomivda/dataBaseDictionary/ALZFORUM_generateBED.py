#! /usr/bin/env python3
import time, argparse #, mariadb
from dataBaseDictionary.ChromosomeToGenePosition_Dictionary import dictionaryFromBED
## -f ALZFORUM_AlzGene.txt -c 1 -s ':' -d /home/flanduzzi@iit.local/DataBases_BioInfo/Homo_sapiens.gene_info.gz

""" FUNCTION: loadData(filename)
Load the data from the input file.
"""
def loadData(filename, separator="\t", column=0):
    results=[]
    with open(filename,"r") as f:
        line=f.readline()
        while line!="" :
            if line[0]!="#" :
                try :
                    splitted=line.strip().split(separator)
                    results.append(splitted[column])
                except Exception as exc:
                    raise exc
                    pass
            line=f.readline()
    return results

""" FUNCTION: createStandardGeneNameSet(listElements,dictionaryGeneInfo, column=0)
Generate a list with the StandardGeneNames and put the excluded in a secondary list.
"""
def createStandardGeneNameSet(listElements,dictionaryGeneInfo=None, column=0):
    univoqueSet=set()
    discarted=set()
    for element in listElements :
        ### Extract the GeneAliases from the list
        val=None
        if isinstance(element,tuple) or isinstance(element, list):
            val=element[column]
        else :
            val=element
        ### IF a ALIASES Dictionary is AVAILABLE load StandardGeneNames else load name
        if ( dictionaryGeneInfo  is not None ) and ( isinstance(dictionaryGeneInfo, dictionaryGeneInfo_Aliases)  ):
            try:
                ### Obtain the StandardGeneName associated with the Aliases
                uniqueName=dictionaryGeneInfo[val]
                print(f"# {val} --> {uniqueName} ")
                univoqueSet.add(uniqueName)
            except :
                ### If 'val' do not correspond to any ALIASES saved in the DISCARTED
                discarted.add(val)
        else:
            univoqueSet.add(val)
    return univoqueSet,discarted

if(__name__=='__main__'):    
    parser = argparse.ArgumentParser(description="Create a BED file from a list of Genes serching the position in the UCSC Genome Browser. Required the library 'mariadb' to connect with the DataBase (pip install mariadb).")
    # awk '{ if(NR>1){ if($1<3){ print $2; }}}' oncokb_biomarker_drug_associations.csv | sort | uniq | awk '{ printf("%s,", $1) ; }'
    parser.add_argument('-i','--input',action='store',type=str,help="List of genes names separeted by comma ','.", required=False, default=None)
    parser.add_argument('-f','--file',action='store',type=str,help="File written in BED format from which load the record of the dictionary.", required=False, default=None)
    parser.add_argument('-s','--column_separator',action='store',type=str,help="The column separator used in the INPUT FILE (String) [DEFAULT: '\t'].", required=False, default='\t')
    parser.add_argument('-c','--column_gene',action='store',type=int,help="The column in the INPUT FILE that contains the GENE NAME (integer) [DEFAULT: 0].", required=False, default=0)
    
    parser.add_argument('-d','--dictionary',action='store',type=str,help="Dictionary for the aliases of the gene names [EXAMPLE: Homo_sapiens.gene_info.gz, DEFAULT: None]", required=False, default=None) 
    # The Standard Dictionary used for ALIASES is: "/home/flanduzzi@iit.local/DataBases_BioInfo/Homo_sapiens.gene_info.gz"
    
    parser.add_argument('-db','--database',action='store',type=str,help="UCSC Genome Browser database where extract the informations (default 'hg38')", required=False, default='hg38')
    parser.add_argument('-o','--output',action='store',type=str,help="Output file where will be stored the gene information in BED format.", required=False, default=None)
    
    parser.add_argument('-sort','--sorted_output',action='store_true', help="Sort the output results, this will slow down the execution [DEFAULT: False].", default=False)
    parser.add_argument('-v','--verbose',action='store_true',help="Verbose [DEFAULT: False].", default=False)
    arg=parser.parse_args()
    
    ### LOAD ALIASES DICTIONARY
    DICTIONARY=arg.dictionary
    from dataBaseDictionary.GeneInfo_Dictionary import dictionaryGeneInfo_Aliases
    if DICTIONARY is not None :
        start_time = time.time()
        DICTIONARY=dictionaryGeneInfo_Aliases(DICTIONARY)
        print("# GeneInfo Alias{:} - Loading time {:.3f} s".format(DICTIONARY,time.time() - start_time))
        
    ### Manage the two possible INPUT merging them into a single list
    geneList=[]
    if arg.file is not None:
        ### INPUT from FILE
        geneList=loadData(arg.file, separator=arg.column_separator, column=arg.column_gene)
    if arg.input is not None:
        ### INPUT from STRING
        for gene in arg.input.strip().split(","):
            geneList.append(gene.strip())
    
    
    ### VERIFY that at least one gene is loaded
    if len(geneList)>0 : 
        ### CREATE the LIST of GENES (no discarded in this step because there is no DICTIONARY)
        geneSet,discarded=createStandardGeneNameSet(geneList, None)
        
        ### GENERATE the BED DICTIONARY - launching the QUERY on UCSC Genome Browser
        myDictionary=dictionaryFromBED()
        DATABASE=arg.database
        start_time=time.time()
        myDictionary.loadFromUCSC(geneSet, dataBaseName=DATABASE, dictionaryAliases=DICTIONARY, verbose=arg.verbose)
        print("# DictionaryBED[{:d}] {:} - Loading time {:.3f} s".format(myDictionary.size(), myDictionary.source_path, time.time() - start_time))
        
        ### VERIFY which GENEs has been FOUND on 'UCSC Genome Browser' and which NOT
        notFound=geneSet.difference(myDictionary.geneSet)
        notFound=list(notFound)
        notFound.sort()
        print("# Genes not found [{:d}]: {:} ".format(len(notFound), ",".join(notFound)))
        
        if DICTIONARY is not None:
            ### REPEAT the QUERY USING the TRANSLATED NAME from the ALIAS dictionary
            geneSet_alias,notAvailableAliases=createStandardGeneNameSet(notFound, DICTIONARY)
            #discarded=discarded.union(notAvailableAliases)
            #if len(discarded)>0:
            #    discardedList=list(discarded)
            #    discardedList.sort()
            
            availableAliases=set(notFound).difference(notAvailableAliases)
            
            if len(geneSet_alias)>0 :
                myDictionary.loadFromUCSC(geneSet_alias, dataBaseName=DATABASE, dictionaryAliases=DICTIONARY, verbose=arg.verbose)
                print("# DictionaryBED[{:d}] {:} - Loading time {:.3f} s".format(myDictionary.size(), myDictionary.source_path, time.time() - start_time))
                
                finalGeneSet=geneSet.difference(availableAliases).union(geneSet_alias)
                ### VERIFY which GENEs has been FOUND on 'UCSC Genome Browser' and which NOT
                notFound=finalGeneSet.difference(myDictionary.geneSet)
                notFound=list(notFound)
                notFound.sort()
                print("# Genes not found [{:d}]: {:} ".format(len(notFound), ",".join(notFound)))
                
#            ### Check also the names not available in the DICTIONARY - this allow the inclusion of mithocondrial DNA
#            geneSet=geneSet.union(discarded)
        
        
        ### PREPARE for the OUTPUT - on file or on screen
        OUTPUT=arg.output
        if OUTPUT is not None :
            outfile=open(OUTPUT, "w")
            outfile.write(myDictionary.toBED_String(sort_output=arg.sorted_output))
            outfile.close()
        else :
            print(myDictionary.toBED_String(sort_output=arg.sorted_output))
        
    else :
        print(arg)
        print("No Gene selected please see the help message using the '-h' option.")
