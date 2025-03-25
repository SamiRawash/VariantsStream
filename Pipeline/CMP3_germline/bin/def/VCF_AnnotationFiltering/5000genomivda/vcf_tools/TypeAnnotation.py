DEFAULT_HEADER_TypeAnnotation="Gene|GeneID|Annotation|AnnotationImpact|Region|FeatureID|Function|AAmodification|HGVS.c|HGVS.p"
class typeAnnotation :
    
    def __init__(self, gene="", geneID="", meaning="", impact="", region="", featureID="", transcript="", aamodification="", HGVSc="", HGVSp="", 
                 gnomad_af=None, rank=None, cDNA=None, cds=None, aa_position=None):
        self.gene=gene
        self.geneID=geneID
        self.annotation=meaning
        self.annotationImpact=impact
        self.region=region
        self.featureID=featureID
        self.transcript=transcript
        self.AAmodification=aamodification
        self.HGVSc=HGVSc     
        self.HGVSp=HGVSp
        self.gnomAD_AF=gnomad_af
        self.rank=rank
        self.cDNA=cDNA
        self.cds=cds
        self.aa_position=aa_position
        self.isCanonical=None
        
    def __str__(self, verbose=False):
        #subtxt=""
        #if verbose :
        #    subtxt="{:} {:}".format(self.impact, self.region) 
        return "{:}|{:}|{:}|{:}|{:}|{:}|{:}|{:}|{:}|{:}".format(self.gene,self.geneID, self.annotation,self.annotationImpact,self.region,self.featureID,self.transcript,self.AAmodification,self.HGVSc,self.HGVSp)
        
class typeAnnotation_LOF :
    
    def __init__(self, gene="", geneID="", number_transcript=-1, percent=-1):
        self.gene=gene
        self.geneID=geneID
        self.numberOfTranscript=number_transcript
        self.percentOfTranscript=percent
        
    def __str__(self, verbose=False):
        return "{:}|{:}|{:d}|{:.3f}".format(self.gene,self.geneID,self.numberOfTranscript,self.percentOfTranscript)

class typeAnnotation_GnomAD:
    
    def __init__(self, date="", gene="", meaning="", region="", AAmodification="", avisnp="", cosmic="", gnomadAll=-1.,gnomad1=-1., exacAll=-1.,exac1=-1.):
        self.date=date
        self.gene=gene
        #self.geneID=geneID
        self.annotation=meaning
        #self.annotationImpact=impact
        self.region=region
        #self.transcript=transcript
        self.AAmodification=AAmodification
        
        self.avisnp=avisnp
        self.cosmic70=cosmic
        
        self.gnomadAll=gnomadAll
        self.gnomadNFE=gnomad1

        self.exacAll=exacAll
        self.exacNFE=exac1
            
    def loadFieldsFromDictionary(self, dictionary):
        self.date=dictionary['ANNOVAR_DATE']
        self.gene=dictionary['GeneDetail.refGene']
        #self.geneID=geneID
        self.annotation=dictionary['ExonicFunc.refGene']
        #self.annotationImpact=impact
        self.region=dictionary['Func.refGene']
        #self.transcript=transcript
        self.AAmodification=dictionary['AAChange.refGene']
        
        self.avisnp=dictionary['avsnp147']
        self.cosmic70=dictionary['cosmic70']
        
        self.gnomadAll=dictionary['gnomAD_genome_ALL']
        self.gnomadNFE=dictionary['gnomAD_genome_NFE']

        self.exacAll=dictionary['ExAC_ALL']
        self.exacNFE=dictionary['ExAC_NFE']

class typeInfo(dict):
    
    def __init__(self, dictionaryInfo):
        super().__init__()
    
        ### Categories used to interpret the data fields
        self.fieldString=set()
        self.fieldInt=set()
        self.fieldFloat=set()
        self.fieldFlag=set()
        
        self.nonStandardItemSet=set()
        
        # For each field load the corresponding Type
        for key,item in list(dictionaryInfo.items()) :
            if item is None :
                self.nonStandardItemSet.add(key)
            elif item[1].lower() == "string" :
                self.fieldString.add(key)
            elif item[1].lower() == "integer" :
                self.fieldInt.add(key)
            elif item[1].lower() == "float" :
                self.fieldFloat.add(key)
            elif item[1].lower() == "flag" :
                self.fieldFlag.add(key)
            else:
                self.nonStandardItemSet.add(key)

    def __mystr(self,s):
        try:
            return str(s)
        except:
            return None
    
    def __myint(self,s):
        try:
            return int(s)
        except:
            return None
    
    def __myfloat(self,s):
        try:
            return float(s)
        except:
            return None
    
    def __myflag(self,s):
        if s is not None:
            return True
        else:
            return False
    
    def __str__(self):
        txt=""
        for key,value in self.items():
            txt+="{}={};".format(key, value)
        return txt
        
    def cleanData(self):
        #self.storedData=typeAnnotation_GnomAD()
        super().clear()
    
    def extractSubset(self, subSet):
        keys=set(self.keys())-subSet
        return { key : self[key] for key in keys }
            
    def parseText(self, key, value):
        if key in self.fieldString :
            self[key]=self.__mystr(value)
        elif  key in self.fieldInt :
            self[key]=self.__myint(value)
        elif  key in self.fieldFloat :
            self[key]=self.__myfloat(value)
        elif key in self.fieldFlag:
            self[key]=self.__myflag(value)
        else :
            self[key]=self.__mystr(value)