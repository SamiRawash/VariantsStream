import time, mariadb
from vcf_tools.CommonTools import templateCompressReader, chromosomeCode, chromosomeCode2sting

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

    def containsVariant(self, variant):
        #print("!!!!!!! TO BE CHECK !!!!!!! - genomicBED.dictionaryFromBED.containsVariant(variant) ")
        if variant.chromosome in self:
            for gene in self[variant.chromosome]:
                if variant.position < gene[0] :
                    return False
                elif variant.position <= gene[1] :
                    return True
        return False
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
    def loadFromUCSC(self, genesNameSet, dataBaseName='hg38', table='wgEncodeGencodeCompV42', dictionaryAliases=None, verbose=True):
        ### UCSC Settings - Host, Port, Tables NAME
        host="genome-mysql.soe.ucsc.edu"
        port=3306
        #tableCanonical='knownCanonical'
        #tableKgXref='kgXref'
        table="ncbiRefSeq"
        
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
        q="SELECT nrs.chrom, MIN(nrs.txStart), MAX(nrs.txEnd), nrs.name2 FROM {:}.{:} nrs where name2 IN ({:}) GROUP BY nrs.chrom, nrs.name2".format(dataBaseName,table,stringSet[:-1])
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
        newDict=GenomeRegionsDictionary()
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

class GenomeRegionsDictionary(dict):
    pass