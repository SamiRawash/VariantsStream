#! /usr/bin/env python3
'''
Created on Jun 18, 2021

@author: flanduzzi
'''  
### Test INPUT Files
# -i C:/test_data/vcf/T_vs_N_filtered_short_cosmic_annotation.vcf 
# -d C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/dataBaseDictionary/oncokb_biomarker_drug_associations.csv
# -a CSQ -o test_outputPASSING.data
# -g C:/Users/flanduzzi/Workspace_Python/5000GenomiVdA/dataBaseDictionary/Homo_sapiens.gene_info.gz,9606

import argparse
import time
#from dataBaseDictionary.oncoKB_Dictionary import dictionaryOncoKB
#from dataBaseDictionary.GeneInfo_Dictionary import dictionaryGeneInfo
from vcf_tools.CommonTools import chiaveVal, chromosomeCode, chromosomeCode2sting, vcf_columnFORMAT_2String, templateCompressReader
from vcf_tools.parserAnnotationVCF import parserAnnotator_CSQ, parserAnnotator_ANN, parserAnnotator_LOF
from vcf_tools.TypePolymorfism import typePolymorfism_Extended
from vcf_tools.TypeAnnotation import typeInfo,DEFAULT_HEADER_TypeAnnotation

#from dataBaseDictionary.Cosmic_Dictionary import setCosmicCodingMuts
HEADER_Format=['DP','AD','GT'] #"DP\tAD\tGT\t"
HEADER_TypeAnnotation=DEFAULT_HEADER_TypeAnnotation.replace("|", "\t")+"\t"
DEFAULT_INFO_FIELDS="AF,CADD_phred,InterVar_automated,CLNSIG,CLNREVSTAT"

class readerVCF_AnnotationFilter(templateCompressReader):
    
    def __init__(self, source, selectedfieldsINFO=None, output="Table.csv"):
        super().__init__()
        self.source_path=source
        self.vcfType="unknown"
        self.passingValue="PASS"
        self.headerColumns=""
            
        ### Internal Variables
        self.FIELDS_INFO=selectedfieldsINFO
        self.TEMPLATE_INFO=""
        self.outfile=open(output, 'w')
        self.DUMP_SAMPLE_NAME=False
        
        self.KEEP_ONLY_PASS=False
        self.ONLY_TRANSCRIPT_WITH_IMPACT=False
        self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT=False
        
        ### Counters
        self.columnFormat_SampleNumber=0
        ### Sample IDs - Used mainly in VCF from multiple sample: ex. Somatic, Trio, ..
        self.sampleID=None
        
        ### Internal Variable - Header 
        self.dictionaryFormat={}
        self.dictionaryInfo={}
        self.dictionaryFilter={}
        
        ### Key Set and Dictionary of the Parser used to interpret the Data
        self.parserAnnotation=None
        self.annotationParserDictionary={}
        self.setAnnotation=set()
        

    def doWithLines(self, line, variantsCounter):
        ### Increase the variants counter
        variantsCounter = variantsCounter + 1
        variant, filtro = self.line2typePolymorfism(line)
        
        if self.KEEP_ONLY_PASS and filtro!='PASS':
            return variantsCounter
        
        if self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT :
            self.keepOnlyCanonicalTranscripts_SnpEff(variant)
        
        flag_lof=True if 'LOF' in variant.info else False
        flag_nmd=True if 'NMD' in variant.info else False
        
        if 'ANN' in variant.info:
            txt = "{:}\t{:d}\t{:}\t{:}\t{:}\t".format(
                chromosomeCode2sting(variant.chromosome, prefix="chr"), 
                variant.position, variant.idVariant, variant.reference, variant.alteration)
            
            ### FIELDs - FORMAT
            txt_format = vcf_columnFORMAT_2String(variant.format, "\t")
            
            ### FIELDs - INFO
            txt_info = ""
            for key in self.FIELDS_INFO:
                if key in variant.info:
                    val = variant.info[key]
                    txt_info += "{}\t".format(val)
                else:
                    txt_info += "None\t"
            
            #if txt_info[-1] == '\t':
            #    txt_info = txt_info[:-1]
                
            variant_info=variant.info['ANN']
            lista_transcripts=self.selectInterestingTranscripts(variant_info) if self.ONLY_TRANSCRIPT_WITH_IMPACT else variant_info
            for annotation in lista_transcripts:
                newVar = annotation.__str__()
                self.outfile.write("{:}{:}{:}{:}{:}\t{:}\n".format(txt, newVar.replace("|", "\t") + "\t", txt_format, txt_info,flag_lof,flag_nmd))
        
        return variantsCounter
    
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
    
    def loadData(self):
        #file=open(self.source_path,"r")
        file=self.openSourceFile(self.source_path)
        variantsCounter=0
        headerCounter=0
        
        ### Check the Filter Values
        line=self.readLine(file)
        while line != "" :
            if line[0] != "#":
                try:
                    variantsCounter = self.doWithLines(line, variantsCounter)       
                except Exception as exc:
                    raise exc
                    print("EXCEPTION: "+str(exc))
                    print("EXCEPTION RAW: "+line)
            else:    
                if line[:9]=="##FORMAT=" :
                    self.splitHeaderFormat(line)
                elif line[:7]=="##INFO=" :
                    self.splitHeaderInfo(line)
                elif line[:9]=="##FILTER=":
                    self.splitHeaderFilter(line)
                elif line[:13]=="##fileformat=":
                    self.vcfType=line[13:-1]
                    if "VCF" not in self.vcfType.upper():
                        print("Not a VCF format")
                        break
                elif line[0]=="#" and line[:2]!="##" :
                    self.doWithHeaderColumn(line)
                headerCounter=headerCounter+1
            
            #line=file.readline()
            line=self.readLine(file)
        
        # CLOSE the input file
        if not file.closed :
            file.close()
            
        return variantsCounter


    def selectInterestingTranscripts(self, list_transcripts, preference_order={'HIGH':1,'MODERATE':2}):
            
        if list_transcripts is None or ( not isinstance(list_transcripts,set) and not isinstance(list_transcripts,list)):
            return []
        elif self.ONLY_TRANSCRIPT_WITH_IMPACT:
            #list_trs=[]
            #for trs in list_transcripts:
            #    if trs.annotationImpact in preference_order.keys() :
            #       list_trs.append(trs)
                    
            dictTranscripts={}
            for trs in list_transcripts:
                idtrs=trs.geneID
                if (idtrs in dictTranscripts) and (trs.annotationImpact in preference_order.keys()):
                    new_impact=trs.annotationImpact
                    old_impact=dictTranscripts[idtrs].annotationImpact
                    if not (old_impact in preference_order):
                        dictTranscripts[idtrs]=trs
                    elif (new_impact in preference_order) and (preference_order[new_impact] < preference_order[old_impact]):
                        #se new_impact<old_impact
                        dictTranscripts[idtrs]=trs
                elif trs.annotationImpact in preference_order.keys() :
                    dictTranscripts[idtrs]=trs
            list_trs=list(dictTranscripts.values())  
            return list_trs
        else:
            dictTranscripts={}
            for trs in list_transcripts:
                idtrs=trs.geneID
                if idtrs in dictTranscripts :
                    new_impact=trs.annotationImpact
                    old_impact=dictTranscripts[idtrs].annotationImpact
                    if not (old_impact in preference_order):
                        dictTranscripts[idtrs]=trs
                    elif (new_impact in preference_order) and (preference_order[new_impact] < preference_order[old_impact]):
                        #se new_impact<old_impact
                        dictTranscripts[idtrs]=trs
                else:
                    dictTranscripts[idtrs]=trs    
            list_trs=list(dictTranscripts.values())
            #list_trs.sort()  
            return list_trs
    
    def composeTableHeader(self, line):
        ### COLLECT Information about the Sample IDs (relevant expecially for the Somatic, Trio ...)
        self.sampleID=line.strip().split('\t')[9:]
        ### ENUMERATE the SAMPLE in the OUTPUT TABLE
        self.columnFormat_SampleNumber=len(self.sampleID)
        if self.DUMP_SAMPLE_NAME :
            for i in range(self.columnFormat_SampleNumber) :
                self.outfile.write("# Sample{:d}: {:}\n".format(i+1, self.sampleID[i]))
        """
		self.columnFormat_SampleNumber = len(line.split('\t')) - 9
        ### COMPOSE the TEMPLATE and HEADER for the INFO fields - TypeAnnotation SUBGOUP
        txtHeader_Format = ""
        for i in range(self.columnFormat_SampleNumber):
            txtHeader_Format += HEADER_TypeAnnotation

        """
        ### COMPOSE the HEADER for VCF Initial Columns
        txtHeader_std="CHROM\tPOS\tID\tREF\tALT\t"

        ### COMPOSE the TEMPLATE and HEADER for the INFO fields
        txtHeader_Info="\t".join(self.FIELDS_INFO)
        
        ### COMPOSE the TEMPLATE and HEADER for the FORMAT fields
        txtHeader_Format=""
        for el in HEADER_Format : 
            if self.columnFormat_SampleNumber == 1:
                txtHeader_Format+=el+"\t"
            else :
                for i in range(self.columnFormat_SampleNumber):
                    txtHeader_Format+=el+"-{:d}\t".format(i+1)

        return "#{:}{:}{:}{:}\t{:}\t{:}\n".format(txtHeader_std, HEADER_TypeAnnotation, txtHeader_Format, txtHeader_Info,'LoF','NMD') 

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
        self.setAnnotation=set(self.annotationParserDictionary.keys())
        
        headertable=self.composeTableHeader(line)
        self.outfile.write(headertable)
                
    def composeTemplateRow(self,fieldsList,dictionary):
        txt_header=""
        txt_formato=""
        dictionaryField2String={"String":"%({:})s", "Float":"%({:})f", "Flag":"%({:})r", "Integer":"%({:})d"}
        for field in fieldsList :
            if field in dictionary :
                txt_header+=field+" \t"
                fieldData=dictionary[field]
                ss=dictionaryField2String[fieldData[1]]
                txt_formato+=ss.format(field)+"\t"
        return txt_header, txt_formato
    
    def splitHeaderInfo(self, stringa):
        splitted=stringa[8:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        formatNumber=splitted[1][7:]
        formatType=splitted[2][5:]
        description=splitted[3][13:-1]
        
        self.dictionaryInfo[formatId]=(formatNumber,formatType,description)
        #return formatId,(formatNumber,formatType,description)

    
    def splitHeaderFormat(self, stringa):
        splitted=stringa[10:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        formatNumber=splitted[1][7:]
        formatType=splitted[2][5:]
        description=splitted[3][13:-1]
        
        self.dictionaryFormat[formatId]=(formatNumber,formatType,description)
        #return formatId,(formatNumber,formatType,description)
    
    
    def splitHeaderFilter(self, stringa):
        splitted=stringa[8:].strip()[:-1].split(",")
        formatId=splitted[0][3:]
        description=splitted[1][13:-1]
        self.dictionaryFilter[formatId]=description

    
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


    def _chekGeneNameFromAnnotation_Filter(self,variant):
        geneNames=set()
        ### For Each Notation Load the Gene Name in the set of Names
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
        ### For Each Notation Load the Gene Name in the set of Names
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
    def parseInfo(self, stringa):
        splitted=stringa.split(";")
        remeaningAnnotation=set(self.setAnnotation)
        
        # Empty the content of the Data - the values are loaded one by one when the key is encountered
        self.info_field_interpreter.cleanData()
        
        for keyVal in splitted:
            x,y=chiaveVal(keyVal)
            if x in remeaningAnnotation :
                # If is a Formatted Field - Load the specific Parser 
                index,parser=self.annotationParserDictionary[x]
                # Store the parsed data in the properties corresponding the Field
                self.info_field_interpreter[x]=parser.parseTxt(y)

                remeaningAnnotation.remove(x)
            else:
                self.info_field_interpreter.parseText(x,y)

        return self.info_field_interpreter #tuple(results) #affectedGenesSet

    def line2typePolymorfism(self, line):
        ### Parsing della stringa variante in oggetto
        x=line.split("\t") #,maxsplit=9)
    
        chromosome = chromosomeCode(x[0])
        
        posizione=int(x[1])
        idVariant=x[2]
        
        riferimento=x[3]
        alterazione=x[4]

        filtro=x[6]
        informazioni=self.parseInfo(x[7])
        formato=self.parseFieldsFormat(x[8:])
        ### Creazione della tuple contenete i dati 
        variant=typePolymorfism_Extended(chromosome, posizione, idVariant, riferimento, alterazione, informazioni, formato)
        variant.txt=line.strip()
        return variant,filtro
    
    def dictionaryInfo_toString(self):
        txt=""
        for el in self.dictionaryInfo.items() :
            txt+="{:<15} {:12} : {:} \n".format(el[0], "[%s]"%el[1][1], el[1][2])
        return txt
    
    def dictionaryFormat_toString(self):
        txt=""
        for el in self.dictionaryFormat.items() :
            txt+="{:<15} {:12} : {:} \n".format(el[0], "[%s]"%el[1][1], el[1][2])
        return txt
        
    def __str__(self):
        return "[{:}={:d},{:}={:d}] {:}".format("PASS",len(self.listDataPassing),"NotPASS", len(self.listDataNonPassing), self.source_path)


    def close(self):
        if not self.outfile.closed :
            self.outfile.close()


################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="ANNOTATED VCF convert to TABLE: \nconsidering a VCF annotated the script extract a subset of fields selected by the user and print them into a tabular format. ")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file (produced by SPNEff or VEP).", required=True, default=None)
    parser.add_argument('-fi','--fieldsinfo',action='store',type=str,help="List of fields to extract from the field 'INFO' of the VCF (the order matter!) \n[DEFAULT: {:}].".format(DEFAULT_INFO_FIELDS), required=False, default=DEFAULT_INFO_FIELDS)
    parser.add_argument('-s','--dump_sample',action='store_true', help="Dump the name of the sample [DEFAULT: False].", default=False)
    parser.add_argument('-only','--only_transcript_with_impact',action='store_true', help="Flag excluding all transcript that do no have an impact HIGH or MODERATE [DEFAULT: False].", default=False)
    parser.add_argument('-canonical','--keep_only_canonical_transcript',action='store_true', help="Flag excluding all non canonical transcript (is canonical the transcript with the larger CDS or for non coding the larger cDNA) [DEFAULT: False].", default=False)
    
    ### PARAMETERs
    parser.add_argument('-o','--output',action='store',type=str,help="Output file path used to store the tabular results .", required=False, default="Table_output.tsv")
    arg=parser.parse_args()

    INPUT=arg.input
    ### Load the field for INFO
    FIELDS_INFO=None
    if arg.fieldsinfo is not None: 
        FIELDS_INFO=arg.fieldsinfo.strip().split(",")
    print("FIELDS INFO: {:}".format(",".join(FIELDS_INFO)))
        
    #OUTPUT=arg.output
    output=arg.output
    
    print("INPUT FILE: {:}".format(INPUT))
    reader=readerVCF_AnnotationFilter(INPUT, selectedfieldsINFO=tuple(FIELDS_INFO), output=output)
    reader.DUMP_SAMPLE_NAME=arg.dump_sample
    reader.ONLY_TRANSCRIPT_WITH_IMPACT=arg.only_transcript_with_impact
    reader.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT=arg.keep_only_canonical_transcript
    
    ### Loading of the VCF Data
    start_time = time.time()
    counterRow=reader.loadData()
    print("VCF [{:d}] - Loading time {:.3f} s".format(counterRow,time.time() - start_time))
    reader.close()
    
    print("Type: %s"%reader.vcfType)
    print("# Dictionary FORMAT")
    print(reader.dictionaryFormat_toString())
    print()
    print("# Dictionary INFO")
    print(reader.dictionaryInfo_toString())
    print()
    