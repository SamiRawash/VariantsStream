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

from vcf_tools.CommonTools import chromosomeCode2sting, vcf_columnFORMAT_2String
from vcf_tools.TypeAnnotation import DEFAULT_HEADER_TypeAnnotation
from vcf_tools.VCF_Annotation2Table import readerVCF_AnnotationFilter
from vcf_tools.pathogenicity_score_prediction_lib import scorepredictionDict

HEADER_Format=['DP','AD','GT'] #"DP\tAD\tGT\t"
HEADER_TypeAnnotation=DEFAULT_HEADER_TypeAnnotation.replace("|", "\t")+"\t"
DEFAULT_INFO_FIELDS="gnomAD_genome_ALL,MQRankSum,SIFT_score,SIFT_pred,FATHMM_score,FATHMM_pred,CLNSIG,CADD_raw,CADD_phred,DANN_score,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred"

class readerVCF_AnnotationFilter_ScorePrediction(readerVCF_AnnotationFilter):
    
    def __init__(self, source, selectedfieldsINFO=None, output="Table.csv"):
        super().__init__(source, selectedfieldsINFO, output)
        self.dictionary_score_predictor=scorepredictionDict
        set_score_pred=self.dictionary_score_predictor.keys()
        self.list_score_pred=list(set_score_pred)
        self.list_score_pred.sort()
        self.FLAG_PATHOGENIC="D"
        self.FLAG_BENIGN="T"
        self.FLAG_UNKNOWN="U"
        self.DUMP_SAMPLE_NAME=False
        #self.ONLY_TRANSCRIPT_WITH_IMPACT=True
        
    def composeTableHeader(self, line):
        self.columnFormat_SampleNumber = len(line.split('\t')) - 9
        ### ENUMERATE the SAMPLE in the OUTPUT TABLE
        if self.DUMP_SAMPLE_NAME :
            self.outfile.write("# ")
            for i,sampleID in enumerate(line.strip().split('\t')[9:]) :
                self.outfile.write("Sample[{:d}]: {:}; ".format(i+1, sampleID))
            self.outfile.write("\n")
        
        ### COMPOSE the HEADER for VCF Initial Columns
        txtHeader_std = "CHROM\tPOS\tREF\tALT\t"
        ### COMPOSE the TEMPLATE and HEADER for the INFO fields
        txtHeader_Info = "\t".join(self.FIELDS_INFO)

        ### COMPOSE the TEMPLATE and HEADER for the FORMAT fields
        txtHeader_Format=""
        for el in HEADER_Format : 
            if self.columnFormat_SampleNumber == 1:
                txtHeader_Format+=el+"\t"
            else :
                for i in range(self.columnFormat_SampleNumber):
                    txtHeader_Format+=el+"-{:d}\t".format(i+1)

        
        ### COMPOSE the HEADER for the PATHOGENIC SCORE PREDICTION
        txt_score_pred="Pathogenic\tBenign\tUncertain\t"
        for sp in self.list_score_pred:
            txt_score_pred+="sp_"+sp+"\t"
        
        return "#{:}{:}{:}{:}{:}{:}{:}{:}\n".format(txtHeader_std, txt_score_pred, 'LoF\t','NMD\t', 'N.Transc.\t',HEADER_TypeAnnotation, txtHeader_Format, txtHeader_Info) 
        

    def scorePrediction_columns(self, variant):
        txt_score = ""
        counter_pathogenic = 0
        counter_benign = 0
        counter_uncertain = 0
        counter_predictor = 0
        for key_score_pred in self.list_score_pred:
            if key_score_pred in variant.info:
                counter_predictor += 1
                sp = scorepredictionDict[key_score_pred]
                res = sp(variant.info[key_score_pred])
                counter_pathogenic += res[0]
                counter_benign += res[1]
                counter_uncertain += res[2]
                if res[0] == 1:
                    txt_score += self.FLAG_PATHOGENIC + "\t"
                elif res[1] == 1:
                    txt_score += self.FLAG_BENIGN + "\t"
                else:
                    txt_score += self.FLAG_UNKNOWN + "\t"
            else:
                txt_score += ".\t"
        
        txt_score = "{:d}\t{:d}\t{:d}\t{}".format(counter_pathogenic, counter_benign, counter_uncertain, txt_score)
        return txt_score

    def doWithLines(self, line, variantsCounter):
        ### Increase the variants counter
        variantsCounter = variantsCounter + 1
        variant, filtro = self.line2typePolymorfism(line)
        if self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT :
            self.keepOnlyCanonicalTranscripts_SnpEff(variant)
        flag_lof=True if 'LOF' in variant.info else False
        flag_nmd=True if 'NMD' in variant.info else False
        if 'ANN' in variant.info:
            ### FIELDs - VCF Standard Columns 
            txt = "{:}\t{:d}\t{:}\t{:}\t".format(
                chromosomeCode2sting(variant.chromosome, prefix="chr"), 
                variant.position, variant.reference, variant.alteration)
            
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
            if txt_info[-1] == '\t':
                txt_info = txt_info[:-1]
            
            ################# Score Predictor
            txt_score = self.scorePrediction_columns(variant)
            ############################################
            variant_info=variant.info['ANN']
            #variant_info.sort(key=lambda trs:(trs.cds[1] if trs.cds is not None else -1, trs.cDNA[1] if trs.cDNA is not None else -1) , reverse=True)
            lista_transcripts=self.selectInterestingTranscripts(variant_info)
            for annotation in lista_transcripts:
                newVar = annotation.__str__()
                self.outfile.write("{:}{:}{}\t{}\t{:d}\t{:}{:}{:}\n".format(txt,txt_score,flag_lof,flag_nmd,len(variant.info['ANN' ]), newVar.replace("|", "\t") + "\t", txt_format, txt_info))
        
        return variantsCounter

    """
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
    """    
        
        
################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="ANNOTATED VCF convert to TABLE: \nconsidering a VCF annotated the script extract a subset of fields selected by the user and print them into a tabular format. ")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file (produced by SPNEff or VEP).", required=True, default=None)
    parser.add_argument('-fi','--fieldsinfo',action='store',type=str,help="List of fields to extract from the field 'INFO' of the VCF (the order matter!) \n[DEFAULT: {:}].".format(DEFAULT_INFO_FIELDS), required=False, default=DEFAULT_INFO_FIELDS)
    #parser.add_argument('-head','--header-only',action='store_true', help="Read only the header [default False].", default=False)
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
    reader=readerVCF_AnnotationFilter_ScorePrediction(INPUT, selectedfieldsINFO=tuple(FIELDS_INFO), output=output)
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
    