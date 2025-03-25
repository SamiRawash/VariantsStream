#! /usr/bin/env python3
'''
Created on Jun 18, 2021

@author: flanduzzi
'''  
from vcf_tools.CommonTools import chiaveVal, chromosomeCode, templateCompressReader
from vcf_tools.TypeAnnotation import typeInfo
from vcf_tools.TypePolymorfism import typePolymorfism_Extended
from vcf_tools.genomicBED import dictionaryFromBED, GenomeRegionsDictionary
from vcf_tools.parserAnnotationVCF import parserAnnotator_CSQ, parserAnnotator_ANN, parserAnnotator_LOF


VARIANT_IMPACT="HIGH,MODERATE"
REFERENCE_POPULATION='gnomAD_genome_ALL'
MAXIMUM_ALLEL_FREQUENCY=0.05


""" class: readerVCF_AnnotationFilter_BaseClass
A basic version of the class that parse the VCF and selecting a list of GENEs through a BED and some restriction on the ANNOTATION.
- __init__(self, source, selectedGenes, outPass): is the constructor, if other properties has to be added to the class overwrite this method (!!! REMEBER to always CALL the super().__init__(...) !!!)
- doWithVariantLine(self,line): is the function that is called for each VARIANT line in the VCF (usually is not required to overwrite this function).
- doWithPassingAnnotatedVariant(variant): is the function called for the PASSING variant in the regions defined by the BED, so is the candidate to be overwrite to personalize the FILTERING
- close(): function that close the output file, overwrite if more output have been created (!!!REMEMBER to call the super().close()!!!)
"""
class readerVCF_AnnotationFilter_BaseClass(templateCompressReader):
    
    def __init__(self, source, selectedGenes=None, outPass="output_passing.vcf"):
        super().__init__()
        self.source_path=source
        self.vcfType="unknown"
        self.passingValue="PASS"
        self.BLOCKING_ERROR=False
        self.INCLUDE_contig=True
        self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT=False
        
        ### Method to verify the Gene NAME
        self.checkGeneNameFromAnnotation=self._chekGeneNameFromAnnotation_NoFilter
        
        ### VERIFY that the SelectedGenes are the correct Type (accepted values are: None, GenomeRegionsDictionary)
        if  selectedGenes is None or len(selectedGenes)==0 :
            self.selectedGenes=None
        elif isinstance(selectedGenes, GenomeRegionsDictionary) :
            self.selectedGenes=selectedGenes
            self.chromosomeSet=set(selectedGenes.keys())
        else :
            raise TypeError("The 'selectedGeneRegions' must be a GenomeRegionsDictionary!")
            
        ### Internal Variables - that will contain the txt     
        self.headerColumns=""
        self.headerTxt=""
        
        ### Internal Variables - OUTPUT file
        self.outfilePass=open(outPass, 'w')
        
        ### Counters
        self.counterPassing=0        
        
        ### Internal Variable - Header 
        self.dictionaryFormat={}
        self.dictionaryInfo={}
        self.dictionaryFilter={}
        self.dictionaryContigs={}
        
        ### Key Set and Dictionary of the Parser used to interpret the Data
        self.parserAnnotation=None
        self.annotationParserDictionary={}
        self.setAnnotation=set()
        
    def loadData(self):
        #file=open(self.source_path,"r")
        file=self.openSourceFile(self.source_path)
        variantsCounter=0
        headerCounter=0
        
        line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                try:
                    ### VARIANT Lines
                    variantsCounter=variantsCounter+1
                    self.doWithVariantLine(line)
                except Exception as exc:
                    ### IF some ERROR occurs print an exeception
                    print("EXCEPTION: "+str(exc))
                    print("EXCEPTION RAW: "+line)
                    ### IF some ERROR occurs and BOLCKING_ERROR raise the Exception
                    if self.BLOCKING_ERROR:
                        raise exc
            else:
                ### HEADER lines
                line=line.rstrip()
                if line[:9]=="##contig=":
                    if self.INCLUDE_contig :
                        self.headerTxt+=line+'\n'
                    self.splitContigs(line)
                elif line[:9]=="##FORMAT=" :
                    self.headerTxt+=line+'\n'
                    self.splitHeaderFormat(line)
                elif line[:7]=="##INFO=" :
                    self.headerTxt+=line+'\n'
                    self.splitHeaderInfo(line)
                elif line[:9]=="##FILTER=":
                    self.headerTxt+=line+'\n'
                    self.splitHeaderFilter(line)
                elif line[:13]=="##fileformat=":
                    self.headerTxt+=line+'\n'
                    self.vcfType=line[13:-1]
                    if "VCF" not in self.vcfType.upper():
                        print("Not a VCF format")
                        break
                elif line[0]=="#" and line[:2]!="##" :
                    self.headerTxt+=line+'\n'
                    self.doWithHeaderColumn(line)
                else :
                    self.headerTxt+=line+'\n'
                headerCounter=headerCounter+1
            
            #line=file.readline()
            line=self.readLine(file)
        
        # CLOSE the input file
        if not file.closed :
            file.close()
            
        return variantsCounter

    """ FUNCTION: 
    Return True if the VARIANT is in the 'selectedGenes' GROUP  
    """
    def vatiantIsInteresting(self,variant):
        if self.selectedGenes is None:
            return True
        if variant.chromosome in self.selectedGenes :
            for gene in self.selectedGenes[variant.chromosome]:
                if variant.position < gene[0] :
                    return False
                elif variant.position <= gene[1] :
                    return True
        return False
    
    def variantIsInteresting_getGene(self,variant):
        if variant.chromosome in self.selectedGenes :
            for gene in self.selectedGenes[variant.chromosome]:
                if variant.position < gene[0] :
                    return None
                elif variant.position <= gene[1] :
                    return gene[2]
        return None

    """ FUNCTION: doWithVariantLine(self,line)
    Method call for variant line of the VCF file.
    """
    def doWithVariantLine(self,line):
        ### Increase the variants counter
        variant,filtro=self.line2typePolymorfism(line)
        ### Check that the variant PASS the quality Check
        if 'PASS' == filtro and self.vatiantIsInteresting(variant):

            ### Extract the data from the fields: INFO, FORMAT
            variant.info=self.parseInfo(variant.info)
            variant.format=self.parseFieldsFormat(variant.format)
            if self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT :
                self.keepOnlyCanonicalTranscripts_SnpEff(variant)
            ### 
            self.doWithPassingAnnotatedVariant(variant)
            
        return            
    
    """ FUNCTION: keepOnlyCanonicalTranscripts_SnpEff(self, variant)
    Function that apply a filter on the variant object filtering out all non canonical transcript (one transcript for each gene will be kept).
    The canonical transcript is decided using the longest CDS, if not present the longest cDNA and if none are present the impact will be used 'HIGH'>'MODERATE'>'MODIFIER'>'LOW'>'.' .
    """
    def keepOnlyCanonicalTranscripts_SnpEff(self, variant):
        if variant.info is not None and 'ANN' in variant.info :
            variant_info=variant.info['ANN']
            variant_info.sort(key=lambda trs:(trs.cds[1] if trs.cds is not None else -1, trs.cDNA[1] if trs.cDNA is not None else ( {"HIGH":-1,"MODERATE":-2,"MODIFIER":-3,"LOW":-4,".":-5}[trs.annotationImpact] )) , reverse=True)
            filtered_list=[]
            gene_set=set()
            for trs in variant_info :
                if trs.gene not in gene_set :
                    filtered_list.append(trs)
                    gene_set.add(trs.gene)
            variant.info['ANN']=filtered_list
    
    """ FUNCTION: doWithPassingAnnotatedVariant(self, variant)
    Method call for the variant that are PASS and are in the interesting group (a.k.a. in the ranges defined in the BED).
    Notice that the 'variant' is an 'typePolymorfism_Extended' object and the fields INFO and FORMAT are already transformed into dictionary.
    """
    def doWithPassingAnnotatedVariant(self, variant):
        ### DUMP the line in the FILE PASS
        self.outfilePass.write(variant.txt + '\n')
        self.counterPassing += 1

    
    def doWithHeaderColumn(self, line):
        self.headerColumns=line.strip()
        self.info_field_interpreter=typeInfo(self.dictionaryInfo)
        
        # Loading of the Annotation
        if 'CSQ' in self.dictionaryInfo:
            self.annotationParserDictionary['CSQ']=(0,parserAnnotator_CSQ(self.dictionaryInfo['CSQ']))
        if 'ANN' in self.dictionaryInfo:
            self.annotationParserDictionary['ANN']=(1,parserAnnotator_ANN(self.dictionaryInfo['ANN']))
        if 'LOF' in self.dictionaryInfo:
            self.annotationParserDictionary['LOF']=(3,parserAnnotator_LOF(self.dictionaryInfo['LOF']))
        if 'NMD' in self.dictionaryInfo:
            self.annotationParserDictionary['NMD']=(4,parserAnnotator_LOF(self.dictionaryInfo['NMD']))
        """
        if 'ALLELE_END' in self.dictionaryInfo:
            self.annotationParserDictionary['ALLELE_END']=(2,parserAnnotator_GnomAD(self.dictionaryInfo))
        
        """
        
        ### Save the Header in each output VCF
        self.writeHeadersIntoOutput()
        
        ### Generate the SET with ALL possible ANNOTATION
        self.setAnnotation=set(self.annotationParserDictionary.keys())
        return
    
    def writeHeadersIntoOutput(self):
        self.outfilePass.write(self.headerTxt)
    
    ######################################################################
    ### Methods for the HEADER PARSING
    ######################################################################
    """ FUNCTION: splitContigs(self,stringa)
    Using the header row INFO load a new element in the 'dictionaryInfo'.
    """
    def splitHeaderInfo(self, stringa):
        splitted=stringa[8:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        formatNumber=splitted[1][7:]
        formatType=splitted[2][5:]
        description=splitted[3][13:-1]
        
        self.dictionaryInfo[formatId]=(formatNumber,formatType,description)
        return
    
    """ FUNCTION: splitHeaderFormat(self,stringa)
    Using the header row FORMAT load a new element in the 'dictionaryFormat'.
    """
    def splitHeaderFormat(self, stringa):
        splitted=stringa[10:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        formatNumber=splitted[1][7:]
        formatType=splitted[2][5:]
        description=splitted[3][13:-1]
        
        self.dictionaryFormat[formatId]=(formatNumber,formatType,description)
        return
    
    """ FUNCTION: splitContigs(self,stringa)
    Using the header row CONTIG load a new element in the 'dictionaryContigs'.
    """
    def splitHeaderFilter(self, stringa):
        splitted=stringa[8:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        description=splitted[1][13:-1]
        self.dictionaryFilter[formatId]=description
        return
    
    """ FUNCTION: splitContigs(self,stringa)
    Using the header row CONTIG load a new element in the 'dictionaryContigs'.
    """
    def splitContigs(self,stringa):
        splitted=stringa.strip()[10:-1].split(",")
        contig=splitted[0].split("=")[1].strip()
        length=int(splitted[1].split("=")[1].strip())
        self.dictionaryContigs[contig]=length
        return
    
    def _chekGeneNameFromAnnotation_Filter(self,variant):
        geneNames=set()
        ### For Each ANNOTATION Load the Gene Name in the set of Names
        for annotation in variant.info['ANN'] :
            if annotation.annotationImpact in self.selectedImpactSet :
                geneNames.add(annotation.gene)
        for annotation in variant.info['CSQ'] :
            if annotation.annotationImpact in self.selectedImpactSet :
                geneNames.add(annotation.gene)
        isOncoKBvariant=bool(geneNames & self.interestingGenesSet)
        return isOncoKBvariant,geneNames

    def _chekGeneNameFromAnnotation_NoFilter(self,variant):
        geneNames=set()
        ### For Each ANNOTATION Load the Gene Name in the set of Names
        for annotation in variant.info['ANN'] :
                geneNames.add(annotation.gene)
        for annotation in variant.info['CSQ'] :
            geneNames.add(annotation.gene)
        if 'Gene.refGene' in variant.info :
            for gene in variant.info['Gene.refGene'].split('\\x') :
                geneNames.add(gene.upper())
        
        isOncoKBvariant=bool(geneNames & self.interestingGenesSet)
        return isOncoKBvariant,geneNames
    
    
    ################################################################
    #### VARIANT ANNOTATION - Parser function
    ################################################################
    """ FUNCTION: parseInfo(self, stringa)
    Populate a dictionary describing the ANNOTATION by extracting the information from the column INFO
    For some special ANNOTATION specified in 'doWithHeaderColumn(self, line)' it use special function to parse them. 
    """
    def parseInfo(self, stringa):
        splitted=stringa.split(";")
        ### To avoid the searching of the annotation already found in the variant it load a set from which remove the name each time they are found.
        remeaningAnnotation=set(self.setAnnotation)
        ### Empty the content of the Data - the values are loaded one by one when the key is encountered
        self.info_field_interpreter.cleanData()
        for keyVal in splitted:
            ### EACH field is managed as KEY,VALUE couple
            x,y=chiaveVal(keyVal)
            
            ### CHECK if the KEY is in the set of SPECIAL ANNOTATION
            if x in remeaningAnnotation :
                ### IF is a SPECIAL ANNOTATION load the corresponding PARSER (also the index is loaded)  
                index,parser=self.annotationParserDictionary[x]
                ### Using the special PARSER store the parsed DATA in the properties corresponding the Field
                self.info_field_interpreter[x]=parser.parseTxt(y)
                ### REMOVE the KEY from the list of SPECIAL ANNOTATION still to be used 
                remeaningAnnotation.remove(x)
            else:
                self.info_field_interpreter.parseText(x,y)

        return self.info_field_interpreter
        
    """ FUNCTION: parseFieldsFormat(self,fieldsFormat)
    Transform the fielts FORMAT into a DICTIONARY where the KEY is the FIELD_NAME and the VALUE is the associated value or tuple (in the case of somatic, trio,... vcf).
    """
    def parseFieldsFormat(self, fieldsFormat):
        tupleSize=len(fieldsFormat)-1
        dictionaryRowVariantFormat={}
        if tupleSize==1 :
            ### Germline case - ONLY ONE COLUMN
            keys=fieldsFormat[0].strip().split(":")
            values=fieldsFormat[1].strip().split(":")
            for n in range(len(keys)):
                dictionaryRowVariantFormat[keys[n].strip()]=(values[n])        
        elif tupleSize==2 :
            ### Normal/Tumor case - TWO COLUMNS
            keys=fieldsFormat[0].strip().split(":")
            values1=fieldsFormat[1].strip().split(":")
            values2=fieldsFormat[2].strip().split(":")
            for n in range(len(keys)):
                dictionaryRowVariantFormat[keys[n].strip()]=(values1[n],values2[n])
        elif tupleSize==0 :
            return None 
        else :
            ### Other case - VARIABLE NUMBER OF FORMAT COLUMNS (slower since uses for loops)
            keys=None
            values=[]
            for n in range(tupleSize+1):
                splittedTxt=fieldsFormat[n].strip().split(":")
                if n==0 :
                    keys=splittedTxt
                else:
                    values.append(splittedTxt)
                    
            for n in range(len(keys)):
                appoggio=[]
                for m in range(tupleSize):
                    appoggio.append(values[m][n])
                dictionaryRowVariantFormat[keys[n].strip()]=tuple(appoggio)

        return dictionaryRowVariantFormat

    """ FUNCTION: line2typePolymorfism(self, line)
    Extract the information from the VARIANT Lines parsing CHROMOSOME, POSITION, REF, VAR, ...
    NOTICE: no parsing of the INFO, FORMAT columns is applied to speedup the execution, the fields are conserved as simple text
    """
    def line2typePolymorfism(self, line):
        ### Parsing della stringa variante in oggetto
        x=line.split("\t")
    
        chromosome = chromosomeCode(x[0])
        
        posizione=int(x[1])
        idVariant=x[2]
        
        riferimento=x[3]
        alterazione=x[4]

        filtro=x[6]
        informazioni=x[7] #self.parseInfo(x[7])
        formato=x[8:] #self.parseFieldsFormat(x[8:])
        
        ### Creazione della tuple contenete i dati 
        variant=typePolymorfism_Extended(chromosome, posizione, idVariant, riferimento, alterazione, informazioni, formato)
        variant.txt=line.strip()
        return variant,filtro
    
    
    ################################################################
    #### METHODS for the DUMPING of the internal DICTIONARY
    ################################################################
    """ FUNCTION: dictionaryInfo_toString(self)
    Transform the dictionary created from the HEADER INFO raws into a STRING to be easily printed.
    """
    def dictionaryInfo_toString(self):
        txt=""
        for el in self.dictionaryInfo.items() :
            txt+="{:<15} {:12} : {:} \n".format(el[0], "[%s]"%el[1][1], el[1][2])
        return txt
    
    """ FUNCTION: dictionaryFORMAT_toString(self)
    Transform the dictionary created from the HEADER FORMAT raws into a STRING to be easily printed.
    """
    def dictionaryFormat_toString(self):
        txt=""
        for el in self.dictionaryFormat.items() :
            txt+="{:<15} {:12} : {:} \n".format(el[0], "[%s]"%el[1][1], el[1][2])
        return txt
    
    """ close(self)
    Close all the output files
    """
    def close(self):
        if not self.outfilePass.closed :
            self.outfilePass.close()
    
    """ FUNCTION: checkImpactInSet(self, variant, selectedImpactSet)
    Check from the SnepEff annotation 'ANN' of the object 'variant' if at least one transcript has IMPACT in the the set 'selectedImpactSet'.
    """
    def checkImpactInSet(self, variant, selectedImpactSet):
        if 'ANN' in variant.info :
            for annotation in variant.info['ANN'] :
                if ( annotation.annotationImpact is not None ) and ( annotation.annotationImpact in selectedImpactSet ):
                    return True
        return False

class readerVCF_AnnotationFilter_InfoFilter(readerVCF_AnnotationFilter_BaseClass):
    
    def __init__(self, source, selectedGenes=None, outPass="output_passing.vcf", filteringCondition=None, filteringPriority=None) : #filteringConditionOnFormat=None, filteringPriority=None, filteringPriorityOnFormat=None):
        super().__init__(source, selectedGenes=selectedGenes, outPass=outPass)
        ### LOAD the DICTIONARY of the FILTERING CONDITION
        if ( filteringCondition is not None ) and ( not isinstance(filteringCondition,dict) ):
            raise TypeError("The 'filteringCondition' must be an instance of 'dict' with key value a function that acts on the corresponding Info field of the VCF.")
        self.filteringCondition=filteringCondition if filteringCondition is not None and len(filteringCondition)>0 else None
        
        ### LOAD the PRIORITY LIST of the FILTERING CONDITION
        self.filteringPriority=None
        if ( filteringPriority is not None ) and ( not isinstance(filteringPriority,list) ):
            raise TypeError("The 'filteringPriority' must be an instance of 'list' with the Info field of the VCF.")
            self.filteringPriority=filteringPriority
        self.filteringPriority=filteringPriority if (( self.filteringCondition is not None) and (filteringPriority is not None )) else list()
        self.FILTERED_FORMAT_FIELDS_IN_INFO=set()
        
    """ OVERWRITE: doWithHeaderColumn(self,line)
    This overwrite manage the 'filteringCondition', 'filteringPriority' synchronization with the VCF 'dictionaryInfo'
    """
    def doWithHeaderColumn(self, line):
        super().doWithHeaderColumn(line)
        if self.filteringCondition is not None :
            ### Remove the FILTER with KEY not present in the VCF DictionaryINFO
            filterKeys=self.filteringCondition.keys() if self.filteringCondition is not None else list()
            vcfInfoKeys=list(self.dictionaryInfo.keys())
            
            for key in self.dictionaryFormat.keys() :
                vcfInfoKeys.append(f"FORMAT.{key}")
                if f"FORMAT.{key}" in self.filteringCondition : 
                    self.FILTERED_FORMAT_FIELDS_IN_INFO.add(f"{key}")
                
            ### if some KEYs of the 'filteringCondition' is not included in the priority
            if len(set(filterKeys)-set(self.filteringPriority)):
                
                ### REMOVE PRIORITY KEYs not present in FILTER CONDITION
                priorityNotPresentInFilter=set(self.filteringPriority).intersection(self.filteringCondition.keys())
                if len(priorityNotPresentInFilter)>0:
                    print("PRIORITIZED KEYS not present in the FILTER CONDITION: ",self._toString_OrderedListFromSet(priorityNotPresentInFilter))
                    print("This keys will be automatically removed!!!")
                    for priorityKey in priorityNotPresentInFilter :
                        self.filteringPriority.remove(priorityKey)
                
                ### ADD MISSING KEYS in the PRIORITY LIST
                missedPriority=(self.filteringCondition.keys())-set(self.filteringPriority)
                if len(missedPriority)>0:
                    missedPriority=list(missedPriority)
                    print("MISSED KEYS in PRIORITY LIST: "," ".join(missedPriority))
                    print("This keys will be automatically added (in alphabetic order) on the tail !!!")
                    for priorityKey in missedPriority :
                        self.filteringPriority.append(priorityKey)
            
                
            ### REMOVE KEYs not present in the VCF
            keyNotPresentInVCF=set(self.filteringPriority)-set(vcfInfoKeys)
            if len(keyNotPresentInVCF)>0:
                print("REMOVE KEY form FILTER CONDITION not present in the VCF: ",self._toString_OrderedListFromSet(keyNotPresentInVCF))
                for priorityKey in keyNotPresentInVCF :
                    self.filteringPriority.remove(priorityKey)
                    self.filteringCondition.pop(priorityKey,None)
            print("PRIORITY LIST: ", " ".join(self.filteringPriority))
            
            
    """ _toString_OrderedListFromSet(self,myset, separator=" ")
    Private function that transform a set in a ordered (alphabetic order) string
    """
    def _toString_OrderedListFromSet(self,myset, separator=" "):
        lista=list(myset)
        lista.sort()
        return separator.join(lista)
    
    """ OVERWRITE: doWithVariantLine(self,line)
    This overwrite manage the filtering through the 'filteringCondition'
    """
    def doWithVariantLine(self,line):
        ### Increase the variants counter
        variant,filtro=self.line2typePolymorfism(line)
        ### Check that the variant PASS the quality Check
        if 'PASS' == filtro and ((self.selectedGenes is None) or self.vatiantIsInteresting(variant)):

            ### Extract the data from the fields: INFO, FORMAT
            variant.info=self.parseInfo(variant.info)
            variant.format=self.parseFieldsFormat(variant.format)
            for key in self.FILTERED_FORMAT_FIELDS_IN_INFO:
                if key in variant.format :
                    value=variant.format[key]
                    filter_key=f"FORMAT.{key}"
                    variant.info[filter_key]=value
            
            variantIsAccepted=True
            ### IF some FILTERING CONDITION has been defined 
            if self.filteringCondition is not None :
                ### Following the PRIORIY list each CONDITION is VERIFIED 
                for key in self.filteringPriority :
                    ### GET the FILTERING FUNCTION
                    filterFunction=self.filteringCondition[key]
                    ### GET the INFO field
                    value=variant.info[key]
                    ### CHECK the CONDITION on the INFO FIELD
                    if not filterFunction(value):
                        variantIsAccepted=False
                        break

            if variantIsAccepted :
                variantIsAccepted=self.crossRequirements(variant)
                
            ### IF ALL CONDITION are VERIFIED the 'doWithPassingAnnotatedVariant(...)' is CALLED   
            if variantIsAccepted :
                self.doWithPassingAnnotatedVariant(variant)
            
        return
    
    def crossRequirements(self, variant) :
        return True
    
class readerVCF_AnnotationFilter(readerVCF_AnnotationFilter_BaseClass):
    
    def __init__(self, source, selectedGenes=None, setImpact=None, 
                 maxAF=MAXIMUM_ALLEL_FREQUENCY, 
                 outPass="output_passing.vcf", 
                 outImpact="output_passing_HighModerate.vcf", 
                 outAF="output_passing_HighModerate_Cosmic_AF.vcf"):
        super().__init__(source, selectedGenes=selectedGenes, outPass=outPass)
        self.field_AF=REFERENCE_POPULATION
        self.maxAF=maxAF
        
        ### Filter on the variation IMPACT
        if setImpact is not None and len(setImpact) :
            self.selectedImpactSet=setImpact
            self.checkGeneNameFromAnnotation=self._chekGeneNameFromAnnotation_Filter
        else:
            self.selectedImpactSet=None
            self.checkGeneNameFromAnnotation=self._chekGeneNameFromAnnotation_NoFilter
        
        ### Internal Variables - OUTPUT file
        self.outfileHighModerate=open(outImpact, 'w')
        self.outfileAllelFrequency=open(outAF, 'w')
        
        ### Counters
        self.counterHighModerate=0
        self.counterAllelFequency=0

    """ Method: doWithPassingAnnotatedVariant(variant)
    This overwrite dump the PASSING variants ('outfilePass') then filter for the IMPACT (dump 'outfileHighModerate') and finally filter for the AF (dump 'outfileAllelFrequency')
    """
    def doWithPassingAnnotatedVariant(self, variant):
        ### DUMP the line in the FILE PASS
        self.outfilePass.write(variant.txt+'\n')
        self.counterPassing+=1
        
        ### VERIFY if the VARIANT have some IMPACTS
        hasIMPACT=False
        if 'ANN' in variant.info :
            for annotation in variant.info['ANN'] :
                if ( annotation.annotationImpact is not None ) and ( annotation.annotationImpact in self.selectedImpactSet ):
                    hasIMPACT=True
                    break
        elif 'CSQ' in variant.info :
            for annotation in variant.info['CSQ'] :
                if ( annotation.annotationImpact is not None ) and ( annotation.annotationImpact in self.selectedImpactSet ):
                    hasIMPACT=True
                    break
        
        ### IF the variant has IMPACT DUMP
        if hasIMPACT :
            
            ### DUMP the line in the FILE HIG-MODERATE
            self.outfileHighModerate.write(variant.txt+'\n')
            self.counterHighModerate+=1
            
            
            ### CHECK the ALLEL FREQUENCY in the desired POPULATION - Notice the Search act also on the CSQ annotation (a.k.a. VEP)
            af=variant.info[self.field_AF] if self.field_AF in variant.info else None
            if af is None and 'CSQ' in variant.info :
                for annotation in variant.info['CSQ'] :
                    af = annotation.gnomAD_AF if af is None or annotation.gnomAD_AF < af else af
                      
            if af is None or af<self.maxAF :
            
                ### IF the FREQUENCY is BELOW maxAF DUMP the line in the FILE AF
                self.outfileAllelFrequency.write(variant.txt+'\n')
                self.counterAllelFequency+=1
    
    def writeHeadersIntoOutput(self):
        super().writeHeadersIntoOutput()
        self.outfileHighModerate.write(self.headerTxt)
        self.outfileAllelFrequency.write(self.headerTxt)

    """ close(self)
    Close all the output files
    """
    def close(self):
        super().close()
        if not self.outfileHighModerate.closed :
            self.outfileHighModerate.close()
        if not self.outfileAllelFrequency.closed :
            self.outfileAllelFrequency.close()


def getGeneticRegions(stringaGeni=None, BEDfile=None, dataBase_hg='hg38', table='wgEncodeGencodeCompV42'):
    newDictionary = dictionaryFromBED()
    ### LOAD the SELECTED GENEs - READINg the BED FILE
    if BEDfile is not None:
        print("BED file to Be LOADED: {:}".format(BEDfile))
        newDictionary.loadBED(BEDfile) ### LOAD the SELECTED GENEs - Reading CHROM,START,STOP from BED
    
    ### LOAD the SELECTED GENEs - UCSC Genome Browser downloading CHROM,START,STOP from BED
    if stringaGeni is not None:
        geneSet = set()
        ### FROM the list of GENEs select ONLY the ONE NOT already LOADED from BED
        for gene in stringaGeni.split(','):
            if not (gene in newDictionary.geneSet):
                geneSet.add(gene) ### IF some GENE from the list is NOT in the dictionary SEARCH on UCSC Genome Browser
        
        if len(geneSet) > 0:
            print("GENEs To Be SEARCHED: {:}".format(", ".join(geneSet)))
            newDictionary.loadFromUCSC(geneSet, dataBaseName=dataBase_hg, table=table)
    print(newDictionary.toBED_String())
    
    return newDictionary.convertSetIntoOrderedLists() if len(newDictionary)>0 else None
    