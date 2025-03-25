import gzip, time

def chiaveVal(txt,separator="="):
    splitted=txt.split(separator, maxsplit=2)
    try:
        return splitted[0],splitted[1]
    except:
        return splitted[0],None

################################### CHROMOSOME 
def chromosomeCode(x, verbose=False):
    if x.startswith("chr"):
        x = x[3:]
    try:
        chromosome = int(x)
    except:
        if x.upper() == 'X':
            chromosome = 100
        elif x.upper() == 'Y':
            chromosome = 101
        elif x.upper() == 'MIT' or x.upper() == 'MT' or x.upper() == 'M':
            chromosome = 102
        else:
            if verbose :
                print("CommonTools NonStandardChromosome:", x)
            chromosome=None
            if x in listNonStandardChromosomes :
                chromosome=listNonStandardChromosomes.index(x)+1000
            else :
                chromosome=len(listNonStandardChromosomes)+1000
                listNonStandardChromosomes.append(x)
    return chromosome

listNonStandardChromosomes=[]
def chromosomeCode2sting(x,prefix=""):
    if x>0 and x<100 :
        return prefix+str(x)
    if x == 100 :
        return prefix+'X'
    elif x == 101 :
        return prefix + 'Y'
    elif x == 102:
        return prefix + 'M'
    elif x>=1000 and x<1000+len(listNonStandardChromosomes):
        return prefix + listNonStandardChromosomes[x-1000]
    return prefix + 'NonStandard'


################################### GENOTYPE
FLAG_HETEROZYGUS="heter"
FLAG_HOMOZYGOUS="homo" 

DICTIONARY_GENOTYPE={
    '0/0':FLAG_HOMOZYGOUS, 
    '0/1':FLAG_HETEROZYGUS, 
    '1/2':FLAG_HETEROZYGUS, 
    '1/1':FLAG_HOMOZYGOUS,
    '0|0':FLAG_HOMOZYGOUS, 
    '0|1':FLAG_HETEROZYGUS, 
    '1|2':FLAG_HETEROZYGUS, 
    '1|1':FLAG_HOMOZYGOUS}

def GenoType2Zygous(genotype):
    splitted=""
    if "/" in genotype :
        splitted=genotype.split("/")
    else :
        splitted=genotype.split("|")
    if len(splitted)==2 :
        return FLAG_HOMOZYGOUS if splitted[0]==splitted[1] else FLAG_HETEROZYGUS 
    else :
        return genotype

def genotype2string(value, separator="\t"):
    txt=""
    if value in DICTIONARY_GENOTYPE :
        txt+=DICTIONARY_GENOTYPE[value]
    else :
        txt+=value
    return txt

def vcf_columnFORMAT_2String(formatColumn, separator="\t"):
    txt=""
    if 'DP' in formatColumn :
        val=formatColumn['DP']
        txt+=separator.join(val) if isinstance(val,tuple) else val
    txt+=separator
    if 'AD' in formatColumn :
        val=formatColumn['AD']
        txt+=separator.join(val) if isinstance(val,tuple) else val
    txt+=separator
    if 'GT' in formatColumn :
        value=formatColumn['GT']
        if isinstance(value,tuple) :
            value=[ genotype2string(x, separator) for x in value ]
            txt+=separator.join(value)
        else :
            txt+=genotype2string(value, separator)
    txt+=separator
    return txt
    
class templateCompressReader:
    
    def __init__(self):
        self.isSourceCompressed=False
    
    def openSourceFile(self, source_path):
        if source_path[-3:] == ".gz":
            file = gzip.open(source_path, 'rb')
            self.isSourceCompressed = True
        else:
            file = open(source_path, "r")
        return file 
       
    def readLine(self,file):
        if self.isSourceCompressed:
            return file.readline().decode()
        else:
            return file.readline()

###### !!!! TO BE TESTED ###################################################################
class templateCompressWriter:
    
    def __init__(self):
        self.isOutputCompressed=False
    
    def openOutputFile(self, output_path):
        if output_path[-3:] == ".gz":
            file = gzip.open(output_path, 'wb')
            self.isSourceCompressed = True
        else:
            file = open(output_path, "w")
        return file 
       
    def writeLine(self,line, file):
        if self.isSourceCompressed:
            return file.write(line.encode())
        else:
            return file.write(line)
