import argparse
import gzip
from vcf_tools.CommonTools import templateCompressReader
"""
Dizionari per tradurre i diversi nomi dei geni in un formato comune
dictionaryGeneInfo: gene -> aliases
dictionaryGeneInfo_Aliases: alias -> gene
"""

class dictionaryGeneInfo(dict, templateCompressReader):
    
    def __init__(self, source, idFilter=None):
        super(dict, self).__init__()
        super(templateCompressReader, self).__init__()
        
        self.source_path=source
        self.idFilter=None if idFilter is None else str(idFilter)
        
        self.isSourceCompressed=False
        file = self.openSourceFile(self.source_path)
        
        line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                splitted=line.split('\t')
                if self.idFilter is None or splitted[0].strip()==self.idFilter :
                    self[splitted[2].strip().upper()]=set(splitted[4].strip().split("|"))
            
            line=self.readLine(file)
        
        file.close()
        self.geneSet=set(self.keys())
    
    def toTable(self, OUTPUT_TYPE="t"):
        listKey=list(self.keys())
        ### Simple SORT
        listKey.sort()
        ### SORT by Genes Name 
        txt=""
        if OUTPUT_TYPE=="h" :
            txt+="<tr> <th>Gene</th> <th> Aliases</th> </tr>" #<th>Cancer Type</th><th>Number Of Drugs</th><th>Drugs</th></tr> \n"
#        elif OUTPUT_TYPE == "j":
#            import json
#            txt+=""
        else :
            txt+="# Gene\tAliases \n"
            
        for key in listKey:
            aliases=self[key]
            if OUTPUT_TYPE == "h":
                txt+="<tr><td>{:}</td> <td>{:d}</td></tr> \n".format(key,", ".join(aliases))
#            elif OUTPUT_TYPE == "j":
#                drugs=self[key]
#                y = json.dumps(list(drugs))
#                txt+='{'+"'level' : {:d} , 'gene' : '{:}' , 'cancerType': '{:}' , 'numDrugs' : {:} , 'drugs' : {:} , ".format(key[0],key[1],key[2],len(drugs), y)+'} \n'
            else:
                txt+="{:}\t{:} \n".format(key,"|".join(aliases))
        if OUTPUT_TYPE=="h" :
            txt="<table>\n"+txt+"</table>\n"
        
        return txt

    """    
    def __openSourceFile(self):
        if self.source_path[-3:] == ".gz":
            file = gzip.open(self.source_path, 'rb')
            self.isSourceCompressed = True
        else:
            file = open(self.source_path, "r")
        return file 
       
    def __readLine(self,file):
        if self.isSourceCompressed:
            return file.readline().decode().strip()
        else:
            return file.readline()
    """
    def __str__(self):
        return "[{:d}] {:}".format(len(self), self.source_path)

class dictionaryGeneInfo_Aliases(dict,templateCompressReader):
    
    def __init__(self, source, idFilter=None):
        super(dict, self).__init__()
        super(templateCompressReader, self).__init__()
        
        self.source_path=source
        if idFilter is None :
            self.idFilter=None
        elif isinstance(idFilter, str) :
            self.idFilter=idFilter
        else :
            self.idFilter=str(idFilter)
        self.geneSet=set()
        
        self.isSourceCompressed=False
        #file = self.__openSourceFile()
        file = self.openSourceFile(self.source_path)
        
        #line=self.__readLine(file)
        line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                splitted=line.split('\t')
                if self.idFilter is None or splitted[0]==self.idFilter :
                    gene=splitted[2].strip().upper()
                    self.geneSet.add(gene)
                    if not splitted[4].strip()=="-" :
                        for alias in splitted[4].strip().split("|") :
                            self[alias.strip().upper()]=gene
                    self[gene]=gene
            #line=self.__readLine(file)
            line=self.readLine(file)
        
        file.close()

    def toTable(self, OUTPUT_TYPE="t"):
        listKey=list(self.keys())
        ### Simple SORT
        listKey.sort()
        ### SORT by Genes Name 
        txt=""
        if OUTPUT_TYPE=="h" :
            txt+="<tr> <th>Alias</th> <th>Gene</th> </tr>" #<th>Cancer Type</th><th>Number Of Drugs</th><th>Drugs</th></tr> \n"
#        elif OUTPUT_TYPE == "j":
#            import json
#            txt+=""
        else :
            txt+="# Alias\tGene \n"
            
        for alias in listKey:
            gene=self[alias]
            if OUTPUT_TYPE == "h":
                txt+="<tr><td>{:}</td> <td>{:d}</td></tr> \n".format(alias,gene)
#            elif OUTPUT_TYPE == "j":
#                drugs=self[key]
#                y = json.dumps(list(drugs))
#                txt+='{'+"'level' : {:d} , 'gene' : '{:}' , 'cancerType': '{:}' , 'numDrugs' : {:} , 'drugs' : {:} , ".format(key[0],key[1],key[2],len(drugs), y)+'} \n'
            else:
                txt+="{:}\t{:} \n".format(alias,gene)
        if OUTPUT_TYPE=="h" :
            txt="<table>\n"+txt+"</table>\n"
        
        return txt

    def __str__(self):
        return "[{:d}] {:}".format(len(self), self.source_path)


if(__name__=='__main__'):
    import time
    parser = argparse.ArgumentParser(description="Gene Info - Dictionary")
    parser.add_argument('-i','--input',action='store',type=str,help="Input GENE INFO data file.", required=False, default="/home/flanduzzi@iit.local/DataBases_BioInfo/Homo_sapiens.gene_info.gz")
    parser.add_argument('-f','--filter',action='store',type=str,help="Filter the column ID by this value (Human=9606)", required=False, default=None)
    parser.add_argument('-l','--list',action='store',type=str,help="List of gene names to be used as example", required=False, default=None)
    #parser.add_argument('-o','--output',action='store',type=str,help="Output VCF data file", required=False, default="C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/actionableGenes/sample2_tumor_normal_mutect2_snpeff.vcf.gz")
    arg=parser.parse_args()
    
    INPUT=arg.input 
    FILTER_ID=arg.filter
    #OUTPUT=arg.output
    
    start_time = time.time()
    GeneInfo_Dictionary=dictionaryGeneInfo(INPUT,FILTER_ID)
    print("GeneInfo{:} - Loading time {:.3f} s".format(GeneInfo_Dictionary,time.time() - start_time))
#    print(myDictionary.toTable())
    
    start_time = time.time()
    GeneInfoAlias_Dictionary=dictionaryGeneInfo_Aliases(INPUT,FILTER_ID)
    print("GeneInfo Alias{:} - Loading time {:.3f} s".format(GeneInfoAlias_Dictionary,time.time() - start_time))
    
    if arg.list is not None :
        GENES=arg.list.upper().split(",")
        GENES=set(GENES)
        
        print("### GENE Aliases")
        for gene in GENES :
            if gene in GeneInfo_Dictionary:
                print("%-10s\t%s"%(gene,", ".join(GeneInfo_Dictionary[gene])))
        
        
    
#    print(myDictionary_alias.toTable())