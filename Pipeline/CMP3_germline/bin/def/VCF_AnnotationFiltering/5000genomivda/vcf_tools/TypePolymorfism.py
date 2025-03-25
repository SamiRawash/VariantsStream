from vcf_tools.CommonTools import chromosomeCode, chromosomeCode2sting

class typePolymorfism :

    def __init__(self, chrom="unknown", pos_start=0, idVariant="unknown", ref="unknown", alt="unknown",info="", formato=None):
        self.chromosome=chrom
        self.position=pos_start
        self.idVariant=idVariant
        self.reference=ref
        self.alteration=alt
        self.info=info
        self.format=formato
    
    def chromosome2string(self):
        return chromosomeCode2sting(self.chromosome)
    """
        if self.chromosome<100:
            return str(self.chromosome)
        if self.chromosome==100:
            return 'X'
        elif self.chromosome==101:
            return 'Y'
        elif self.chromosome==102:
            return 'M'
        else:
            return 'unknown'
    """
    def parseLineVCF(self,line):
        ### Parsing della stringa variante in oggetto
        x=line.split("\t") #,maxsplit=9)
        self.chromosome = chromosomeCode(x[0])
        self.position=int(x[1])
        self.reference=x[3]
        self.alteration=x[4]
        self.info=x[7]
        
    def __str__(self):
        return "{:}\t{:d}\t{:}\t{:}\t{:}".format(self.chromosome2string(),self.position,self.reference,self.alteration,self.info)


class typePolymorfism_Extended(typePolymorfism):
    
    def __init__(self, chrom="unknown", pos_start=0, idVariant="unknown", ref="unknown", alt="unknown",info="", formato=None, txtRow=None):
        super().__init__(chrom, pos_start, idVariant, ref, alt, info, formato)
        self.txt=None if txtRow is None else txtRow.strip()
    
    def __str__(self):
        if isinstance(self.info,dict):
            subtxt=self.info.__str__() 
        return "{:}\t{:d}\t{:}\t{:}\t{:}\t{:}".format(self.chromosome2string(),self.position,self.idVariant,self.reference,self.alteration,subtxt)
