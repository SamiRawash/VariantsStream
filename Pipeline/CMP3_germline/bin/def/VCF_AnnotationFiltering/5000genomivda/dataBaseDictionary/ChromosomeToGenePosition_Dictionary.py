import time, argparse, mariadb
from vcf_tools.CommonTools import templateCompressReader, chromosomeCode, chromosomeCode2sting
"""
To install MariaDB the ConnectorC is required: 
curl -sS https://downloads.mariadb.com/MariaDB/mariadb_repo_setup | sudo bash
sudo apt -y install libmariadb3 libmariadb-dev
pip3 install mariadb
"""
class dictionaryFromBED(dict, templateCompressReader):

    def __init__(self):
        super(dict, self).__init__()
        super(templateCompressReader, self).__init__()
        self.source_path=None
        self.isSourceCompressed=False
        self.geneSet=set()
    
    
    """
    dictionary.loadBED : function that load a BED file into the dictionary
    
    UCSC Genome Browser - BED File Format
    The first three required BED fields are:

        1)  chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
        2)  chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
        3)  chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature, however, the number in position format will be represented. For example, the first 100 bases of chromosome 1 are defined as chrom=1, chromStart=0, chromEnd=100, and span the bases numbered 0-99 in our software (not 0-100), but will represent the position notation chr1:1-100. Read more here.
        
    Additional optional BED fields are:
        
        4)  name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode.
        5)  score - A score between 0 and 1000. If the track line useScore attribute is set to 1 for this annotation data set, the score value will determine the level of gray in which this feature is displayed (higher numbers = darker gray). This table shows the Genome Browser's translation of BED score values into shades of gray:
        6)  shade score in range      ≤ 166    167-277    278-388    389-499    500-611    612-722    723-833    834-944    ≥ 945
        7)  strand - Defines the strand. Either "." (=no strand) or "+" or "-".
        8)  thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). When there is no thick part, thickStart and thickEnd are usually set to the chromStart position.
        9)  thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
        10) itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). If the track line itemRgb attribute is set to "On", this RBG value will determine the display color of the data contained in this BED line. NOTE: It is recommended that a simple color scheme (eight colors or less) be used with this attribute to avoid overwhelming the color resources of the Genome Browser and your Internet browser.
        11) blockCount - The number of blocks (exons) in the BED line.
        12) blockSizes - A comma-separated list of the block sizes. The number of items in this list should correspond to blockCount.
        13) blockStarts - A comma-separated list of block starts. All of the blockStart positions should be calculated relative to chromStart. The number of items in this list should correspond to blockCount.
        
    """

    """
    LOAD Data from BED file "chromosome stratPos endPos GeneName ..."
    """
    def loadBED(self, source, fieldSeparator='\t', discardFirstRow=False):
        self._addNameToSources(source)
        counterAddedLines=0
        ### Open the File and Check if is Compressed
        self.isSourceCompressed=False
        file = self.openSourceFile(source)
        
        ### Reading the SOURCE and eventually discard the Header Row
        line=self.readLine(file)
        if discardFirstRow :
            line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                x=line.strip().split(fieldSeparator)
                try :
                    self._addElement(x)
                except:
                    print("EXCEPTION RowNotLoaded: {:}".format(line.strip()))
                counterAddedLines+=1
            line=self.readLine(file)
        
        file.close()
        return counterAddedLines

    """
    LOAD Data from a list of Tuples similar to BED format (chromosome, stratPos, endPos, GeneName, ...)
    """
    def loadList(self, lista):
        self._addNameToSources("List")
        counterAddedLines=0
        for x in lista:
            try :
                self._addElement(x)
            except:
                print("EXCEPTION ElementNotLoaded: ",x)
            counterAddedLines+=1
        return counterAddedLines
    
    """
    LOAD Data from the UCSC Genome Browser:
    it allows also the interpretation of the names through a dictionary of aliases ( alias -> geneSymbol )
    """
    def loadFromUCSC(self, genesNameSet, dataBaseName='hg38', dictionaryAliases=None, verbose=True):
        ### UCSC Settings - Host, Port, Tables NAME
        host="genome-mysql.soe.ucsc.edu"
        port=3306
        #tableCanonical='knownCanonical'
        #tableKgXref='kgXref'
        tableNcbiRefSeq='wgEncodeGencodeCompV42'
        
        ###
        self._addNameToSources(host+":"+str(port)+":"+dataBaseName)
        counterAddedLines=0
        
        # CREATE the string with genes name used in the query
        stringSet=""
        for name in genesNameSet :
            stringSet+="'{:}',".format(name)
        
        start_time = time.time()
        ucscGenomeBrowser=mariadb.connect(host=host, user="genome", port=port, database=dataBaseName)
        cursor=ucscGenomeBrowser.cursor()
        if verbose:
            print("UCSC Genome Browser - Establishing connection time {:.3f} s".format(time.time() - start_time))
        
        #q="SELECT c.chrom, c.chromStart, c.chromEnd, g.geneSymbol from {:}.{:} g inner join {:}.{:} c on g.kgID = c.transcript WHERE geneSymbol IN ({:})".format(dataBaseName, tableKgXref, dataBaseName, tableCanonical, stringSet[:-1])
        q="SELECT nrs.chrom, MIN(nrs.txStart), MAX(nrs.txEnd), nrs.name2 FROM {:}.{:} nrs where name2 IN ({:}) GROUP BY nrs.chrom, nrs.name2".format(dataBaseName,tableNcbiRefSeq,stringSet[:-1])
        if verbose :
            print("execute: "+q)
        cursor.execute(q)
        for row in cursor:
            try :
                self._addElement(row)
            except:
                print("EXCEPTION RowNotLoaded: ",row)
            counterAddedLines+=1
        
        return counterAddedLines

    ### PROTECTED FUNCTION - Required to LOAD Data
    def _addNameToSources(self, newSource):
        if self.source_path is None :
            self.source_path=newSource
        else:
            self.source_path+="|"+newSource
            
    def _addElement(self, x):
        chrom = chromosomeCode(x[0])
        geneName=x[3] if len(x) > 3 and x[3] is not None else None
        if chrom in self:
            self[chrom].add((int(x[1]), int(x[2]), geneName))
        else:
            newSet = set()
            newSet.add((int(x[1]), int(x[2]), geneName))
            self[chrom] = newSet
        if geneName is not None :
            self.geneSet.add(geneName)
            
    
    """def toTable(self, OUTPUT_TYPE="t"):
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
    """ Function: size
    Count the total number of records contained in the Dictionary (included the item in each subset)
    """
    def size(self):
        counter=0
        for key,val in self.items() :
            counter+=len(val)
        return counter

    """ Function: convertSetIntoOrderedLists
    Transform the object into a new dictionary that have as elements the SORTED LIST of gene ranges
    It is useful for a faster search of a pointlike modification.
    """
    def convertSetIntoOrderedLists(self):
        newDict={}
        for key,val in self.items() :
            ordList=list(val)
            ordList.sort()
            newDict[key]=ordList
        return newDict

    """ Function: toBED_String
    Transform the object into a text formatted as BED file.
    """
    def toBED_String(self, sort_output=False):
        txt="#chrom\tchromStart\tchromEnd\tname\n"
        chrList=list(self.keys())
        chrList.sort()
        for x in chrList :
            element_list=sorted(self[x], key=lambda tup: (tup[0],tup[1]) ) if sort_output else self[x]
            for rangeVal in element_list :
                txt+="{:}\t{:d}\t{:d}{:}".format(chromosomeCode2sting(x, 'chr'),rangeVal[0],rangeVal[1],"\n" if rangeVal[2] is None else "\t{:}\n".format(rangeVal[2]) )
        return txt
    
    def __str__(self):
        return "[{:d}] {:}".format(len(self), self.source_path)

class dictionaryCancerGeneCensus(dict, templateCompressReader):

    def __init__(self):
        super(dict, self).__init__()
        super(templateCompressReader, self).__init__()
        self.source_path=None
        self.isSourceCompressed=False
        self.geneSet=set()
    
    """
    LOAD Data from BED file "chromosome stratPos endPos GeneName ..."
    """
    def loadCSV(self, source, fieldSeparator='\t', discardFirstRow=False):
        self._addNameToSources(source)
        counterAddedLines=0
        ### Open the File and Check if is Compressed
        self.isSourceCompressed=False
        file = self.openSourceFile(source)
        
        ### Reading the SOURCE and eventually discard the Header Row
        line=self.readLine(file)
        if discardFirstRow :
            line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                x=line.strip().split(fieldSeparator)
                
                newList=[]
                counter=0
                openTextField=False
                
                for field in x :
                    
                    if not openTextField :
                        counter+=1
                        newList.append(field)
                    else :
                        newList[counter-1]+=field

                    if field.count('"')%2 != 0 :
                        openTextField = not openTextField
                    
                try :
                    self._addElement(newList)
                except:
                    print("EXCEPTION RowNotLoaded: {:}".format(line.strip()))
                    
                counterAddedLines+=1
            line=self.readLine(file)
        
        file.close()
        return counterAddedLines

    ### PROTECTED FUNCTION - Required to LOAD Data
    def _addNameToSources(self, newSource):
        if self.source_path is None :
            self.source_path=newSource
        else:
            self.source_path+="|"+newSource
            
    def _addElement(self, x):
#Gene Symbol,Name,Entrez GeneId,Genome Location,Tier,Hallmark,Chr Band,Somatic,Germline,Tumour Types(Somatic),Tumour Types(Germline),Cancer Syndrome,Tissue Type,Molecular Genetics,Role in Cancer,Mutation Types,Translocation Partner,Other Germline Mut,Other Syndrome,Synonyms
#A1CF,APOBEC1 complementation factor,29974,10:50799421-50885675,2,,11.23,yes,,melanoma,,,E,,oncogene,Mis,,,,"29974,A1CF,ACF,ACF64,ACF65,APOBEC1CF,ASP,ENSG00000148584.14"
        #try :
        posInfo=x[3].strip().split(':')
        chrom = chromosomeCode(posInfo[0])
        posInfo=posInfo[1].strip().split('-')
        startPos=int(posInfo[0])
        endPos=int(posInfo[1])
        geneName=x[0]
        if chrom in self:
            self[chrom].add((startPos, endPos, geneName))
        else:
            newSet = set()
            newSet.add((startPos, endPos, geneName))
            self[chrom] = newSet
        if geneName is not None :
            self.geneSet.add(geneName)
        #except:
        #    raise Exception("Error parsing row")

    """ Function: size
    Count the total number of records contained in the Dictionary (included the item in each subset)
    """
    def size(self):
        counter=0
        for val in self.values() :
            counter+=len(val)
        return counter

    """ Function: convertSetIntoOrderedLists
    Transform the object into a new dictionary that have as elements the SORTED LIST of gene ranges
    It is useful for a faster search of a pointlike modification.
    """
    def convertSetIntoOrderedLists(self):
        newDict={}
        for key,val in self.items() :
            ordList=list(val)
            ordList.sort()
            newDict[key]=ordList
        return newDict

    """ Function: toBED_String
    Transform the object into a text formatted as BED file.
    """
    def toBED_String(self):
        txt="#chrom\tchromStart\tchromEnd\tname\n"
        chrList=list(self.keys())
        chrList.sort()
        for x in chrList :
            for rangeVal in self[x] :
                txt+="{:}\t{:d}\t{:d}{:}".format(chromosomeCode2sting(x, 'chr'),rangeVal[0],rangeVal[1],"\n" if rangeVal[2] is None else "\t{:}\n".format(rangeVal[2]) )
        return txt
    
    def __str__(self):
        return "[{:d}] {:}".format(len(self), self.source_path)


if(__name__=='__main__'):    
    from dataBaseDictionary.GeneInfo_Dictionary import dictionaryGeneInfo_Aliases
    parser = argparse.ArgumentParser(description="Create the Dictionary - UCSC Genome Browser")
    # awk '{ if(NR>1){ if($1<3){ print $2; }}}' oncokb_biomarker_drug_associations.csv | sort | uniq | awk '{ printf("%s,", $1) ; }'
    parser.add_argument('-i','--input',action='store',type=str,help="List of genes names separeted by comma ','.", required=False, default=None)
    parser.add_argument('-f','--file',action='store',type=str,help="File written in BED format from which load the record of the dictionary.", required=False, default=None)
    parser.add_argument('-d','--dictionary',action='store',type=str,help="Dictionary for the aliases of the gene names", required=False, default="C:/Users/flanduzzi/Bioinformatic_DataBase/Homo_sapiens.gene_info.gz") #"C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/dataBaseDictionary/Homo_sapiens.gene_info.gz")
    parser.add_argument('-db','--database',action='store',type=str,help="UCSC Genome Browser database where extract the informations (default 'hg38')", required=False, default='hg38')
    parser.add_argument('-o','--output',action='store',type=str,help="Output file where will be stored the gene information in BED format.", required=False, default=None)
    
    arg=parser.parse_args()
    


    ### LOAD ALIASES DICTIONARY
    DATABASE=arg.database
    DICTIONARY=arg.dictionary
    start_time = time.time()
    DICTIONARY=dictionaryGeneInfo_Aliases(DICTIONARY)
    print("GeneInfo Alias{:} - Loading time {:.3f} s".format(DICTIONARY,time.time() - start_time))
    
    print("######### Section - UCSC Genome Browser")
    ### INPUT from UCSC Genome Browser - Selected from a Gene Symbol List
    geneSet=arg.input
    if geneSet is not None and geneSet.replace(",","").strip()!="" :
        lista=geneSet.strip().replace(" ","").split(',')
        geneSet=set()
        for alias in lista:
            if DICTIONARY is None :
                geneSet.add(alias.strip())
            elif alias.strip() in DICTIONARY :
                geneSet.add(DICTIONARY[alias.strip()])
    else :
        geneSet=None
    
    ### INPUT from a BED File
    sourceFile=arg.file
    
    ### CHECK : at least one input must be specified
    if sourceFile is None and geneSet is None : 
        raise Exception("No imput are set the program require at least one either a source BED or genes list to be searched in UCSC Genome Browser.")
    
    start_time = time.time()
    myDictionary=dictionaryFromBED()
    ### Loading from UCSC
    if geneSet is not None :
        myDictionary.loadFromUCSC(geneSet, dataBaseName=DATABASE, dictionaryAliases=DICTIONARY, verbose=True)
    ### Loading from BED
    if sourceFile is not None :
        myDictionary.loadBED(sourceFile, '\t', True)
    
    print("DictionaryBED[{:d}] {:} - Loading time {:.3f} s".format(myDictionary.size(), myDictionary.source_path, time.time() - start_time))
    sortedDict=myDictionary.convertSetIntoOrderedLists()
    chrList=list(sortedDict.keys())
    chrList.sort()
    for x in chrList :
        print("Chromosome {}".format(chromosomeCode2sting(x)))
        for rangeVal in sortedDict[x] :
            print("\t - {:}[{:d}] : [{:d}:{:d}]".format(rangeVal[2],rangeVal[1]-rangeVal[0], rangeVal[0],rangeVal[1]))
            
    ### Creation of the BUFFER for the OUTPUT file
    OUTPUT=arg.output
    if OUTPUT is not None :
        outfile=open(OUTPUT, "w")
        outfile.write(myDictionary.toBED_String())
        outfile.close()
    
    print()
    print("######### Section - CancerGeneCensus")
    cgcDictionary=dictionaryCancerGeneCensus()
    start_time = time.time()
    cgcDictionary.loadCSV("C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/dataBaseDictionary/cancer_gene_census.csv", fieldSeparator=',', discardFirstRow=True)
    print("DictionaryCancerGeneCensus[{:d}] {:} - Loading time {:.3f} s".format(cgcDictionary.size(), cgcDictionary.source_path, time.time() - start_time))
    ### Creation of the BUFFER for the OUTPUT file
    if OUTPUT is not None :
        outfile=open("C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/dataBaseDictionary/cancer_gene_census.bed", "w")
        outfile.write(cgcDictionary.toBED_String())
        outfile.close()
    
    
    
    print()
    print("######### Section - Intersection")
    setCGC=set()
    for gene in cgcDictionary.geneSet :
        if gene in DICTIONARY :
            setCGC.add(DICTIONARY[gene])
        else :
            print("CGC  NotFound: ",gene)
    print("CancerGeneCensus: {:d}".format(len(setCGC)))    
    
    setUCSC=set()
    for gene in myDictionary.geneSet :
        if gene in DICTIONARY :
            setUCSC.add(DICTIONARY[gene])
        else :
            print("UCSC NotFound: ",gene)
    print("UCSCGenomeBrowser: {:d}".format(len(setUCSC)))    
    
    
    setIntersection=setCGC.intersection(setUCSC)
    listaOrd=list(setIntersection)
    listaOrd.sort()
    print("Intersection of the two sets [{:d}]: \n{:}".format(len(setIntersection), ",".join(listaOrd)))

    setOncoKBNotCGC=setUCSC-setIntersection
    listaOrd=list(setOncoKBNotCGC)
    listaOrd.sort()
    print("OncoKB not in CancerGeneCensus [{:d}]: \n{:}".format(len(setOncoKBNotCGC), ",".join(listaOrd)))
    
    
    
