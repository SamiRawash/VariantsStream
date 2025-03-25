#! /usr/bin/env python3
'''
Created on Jun 18, 2021

@author: flanduzzi
'''  
import argparse, time
from vcf_tools.VCF_AnnotationFilter import readerVCF_AnnotationFilter_InfoFilter, getGeneticRegions


MINIMUM_DEPTH=10
MINIMUM_AF=0.2
VARIANT_IMPACT="HIGH,MODERATE"
REFERENCE_POPULATION='gnomAD_genome_ALL'
MAXIMUM_GNOMAD=0.01

MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS=4
CLINVAR_benign=set(["benign","benign/likely_benign","likely_benign"])
CLINVAR_uncertain=set(["uncertain_significance", ".", "conflicting_interpretations_of_pathogenicity", "not_provided", "conflicting_data_from_submitters", "other"])
FATHMM_maximum=-1.5
SIFT_maximum=0.05
POLYPHEN2_HDIV_minimum=0.453 #0.908 on the FOREUM project
POLYPHEN2_HVAR_mimimum=0.447 #0.908 on the FOREUM project
CADD_minimum=20
DANN_minimum=0.90
GERP_minimum=2


def score_predictor_FATHMM_score(value):
    if value is None or value==".":
        return 0
    try: 
        if float(value)<=FATHMM_maximum:
            return 1
        else:
            return  0
    except ValueError:
        return 0
    return 0
    
def score_predictor_SIFT_score(value):
    if value is None or value==".":
        return 0
    try: 
        if float(value)<=SIFT_maximum:
            return 1
        else:
            return  0
    except ValueError:
        return 0
    return 0
    
def score_predictor_POLYPHEN2_HDIV_score(value):
    if value is None or value==".":
        return 0
    try: 
        if float(value)>=POLYPHEN2_HDIV_minimum:
            return 1
        else:
            return  0
    except ValueError:
        return 0
    return 0

def score_predictor_POLYPHEN2_HVAR_score(value):
    if value is None or value==".":
        return 0
    try:
        if float(value)>=POLYPHEN2_HVAR_mimimum:
            return 1
        else:
            return 0
    except ValueError:
        return 0
    return 0

def score_predictor_CADD_score(value):
    if value is None or value==".":
        return 0
    try:
        if float(value)>=CADD_minimum:
            return 1
        else:
            return 0
    except ValueError:
        return 0
    return 0

def score_predictor_DANN_score(value):
    if value is None or value==".":
        return 0
    try: 
        if float(value)>=DANN_minimum:
            return 1
        else:
            return 0
    except ValueError :
        return 0
    return 0

def score_predictor_GERP(value):
    if value is None or value==".":
        return 0
    try:
        if float(value)>=GERP_minimum:
            return 1
        else:
            return 0
    except ValueError:
        return 0
    return 0

def score_predictor_ClinPred(value):
    if value is not None and value=="D":
        return 1
    return 0

def score_predictor_PROVEAN(value):
    if value is not None and value=="D":
        return 1
    return 0


class readerVCF_AnnotationFilter_Report(readerVCF_AnnotationFilter_InfoFilter):
    
    def __init__(self, source, selectedGenes=None, outPass="output_passing.vcf", filteringCondition=None, filteringPriority=None):
        super().__init__(source, selectedGenes=selectedGenes, outPass=outPass, filteringCondition=filteringCondition, filteringPriority=filteringPriority)
        
        ### Filter on the variation IMPACT
        self.selectedImpactSet=set(VARIANT_IMPACT.split(","))
        self.MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS=MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS
        
    def crossRequirements(self,variant):

        ### CHECK the IMPACT on SNPEFF and VEP annotation
        snpeff=variant.info['ANN'] if 'ANN' in variant.info else []
        vep=variant.info['CSQ'] if 'CSQ' in variant.info else []
        hasIMPACT,variant.info['ANN'],variant.info['CSQ']=self.hasIMPACT(snpeff,vep)
        if not hasIMPACT:
            return False
        
        ### CHECK INHERITANCE
        #clndn=variant.info['CLNDN']
        #gt=variant.format['GT']
        #inheritance=self.doInheritance(clndn,gt)
        #if not inheritance:
        #    return False
        
        ### PATHOGENIC SCORE 
        clinvar=variant.info['CLNSIG']
        if clinvar.lower() in CLINVAR_benign :
            return False
        elif clinvar.lower() in CLINVAR_uncertain :
            FATHMM_score=variant.info['FATHMM_score']
            SIFT_score=variant.info['SIFT_score']
            POLYPHEN2_HDIV_score=variant.info['Polyphen2_HDIV_score']
            POLYPHEN2_HVAR_score=variant.info['Polyphen2_HVAR_score']
            CADD_score=variant.info['CADD_phred']
            DANN_score=variant.info['DANN_score']
            GERP_score=variant.info['GERP++_RS']
            #GERP_nr_score=variant.info['GERP++_NR']
            #GERP_rn_score=variant.info['GERP++_RS_rankscore']
            PROVEAN_pred=variant.info['PROVEAN_pred']
            ClinPred_score=variant.info['ClinPred_pred']
            
            if score_predictor_FATHMM_score(FATHMM_score)+score_predictor_SIFT_score(SIFT_score)+score_predictor_POLYPHEN2_HDIV_score(POLYPHEN2_HDIV_score)+score_predictor_POLYPHEN2_HVAR_score(POLYPHEN2_HVAR_score)+score_predictor_CADD_score(CADD_score)+score_predictor_DANN_score(DANN_score)+score_predictor_GERP(GERP_score)+score_predictor_ClinPred(ClinPred_score)+score_predictor_PROVEAN(PROVEAN_pred)<self.MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS :
                return False
        
        return True
    
    
    def hasIMPACT(self,snpeff,vep):
        hasIMPACT=False
        
        newANN=[]
        for annotation in snpeff :
            if ( annotation.annotationImpact is not None ) and ( annotation.annotationImpact in self.selectedImpactSet ):
                newANN.append(annotation)
                hasIMPACT=True

        newCSQ=[]
        for annotation in vep :
            if ( annotation.annotationImpact is not None ) and ( annotation.annotationImpact in self.selectedImpactSet ):
                newCSQ.append(annotation)
                hasIMPACT=True
        
        return hasIMPACT,newANN,newCSQ
    
    
    ###Function to ADD
    def doInheritance(self, Clinvar_value, GT):
        isValid = False
        if "dominant" in Clinvar_value.casefold() and (GT == "0/1" or GT == "1/0" or GT == "0|1" or GT == "1|0"):
            isValid = True
        elif "recessive" in Clinvar_value.casefold() and (GT == "1/1" or GT == "1|1"):
            isValid = True
        #elif "x-linked" in Clinvar_value.casefold():
        #    gne
        #elif hetcomp??
        return isValid


################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="FILTER ANNOTATED VCF fro REPORT: considering a VCF annotated  with PASS, SNPEFF, ANNOVAR the script select the passing variants that refer to a preselected GENES_SET.")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file (produced by SPNEff, ANNOVAR or VEP).", required=True, default=None)
    parser.add_argument('-blocking','--set_blocking_errors',action='store_true', help="Flag forcing the end of the execution with the raise of an exception in case of an error during the VCF readings [DEFAULT: False].", default=False)
    parser.add_argument('-contigs','--set_include_contigs',action='store_true', help="Flag allowing the inclusion of the CONTIG raws in the VCF header [DEFAULT: False].", default=False)
    
    ### PARAMETERs
    parser.add_argument('-g','--genelist',action='store',type=str,help="List of genes to be searched separeted by commas (no spaces).", required=False, default=None)
    parser.add_argument('-bed','--bedgenefile',action='store',type=str,help="BED file containing the list of genes to be searched.", required=False, default=None)
    parser.add_argument('-impact','--filterimpact',action='store',type=str,help="String containing the variant gene impact separated by commas [DEFAULT: "+VARIANT_IMPACT+"].", required=False, default=VARIANT_IMPACT)
    parser.add_argument('-p','--population',action='store',type=str,help="String describing the field from which extract REFERENCE POPULATION for the allel frequency (AF). Notice that if the information is stored in the 'CSQ' field (a.k.a. VEP annotation) only the population ALL is checked. [DEFAULT: "+REFERENCE_POPULATION+"].", required=False, default=REFERENCE_POPULATION)
    parser.add_argument('-gnomad','--gnomad',action='store',type=float,help="Set the maximum ALLEL FREQUENCY in the reference POPULATION above which discard the variant [DEFAULT: "+str(MAXIMUM_GNOMAD)+"].", required=False, default=MAXIMUM_GNOMAD)
    parser.add_argument('-af','--minAF',action='store',type=float,help="Set the maximum ALLEL FREQUENCY for the ALTERED allele below which the variant is discarded [DEFAULT: "+str(MAXIMUM_GNOMAD)+"].", required=False, default=MAXIMUM_GNOMAD)
    parser.add_argument('-dp','--minDP',action='store',type=float,help="Set minimum DEPTH below which the variant is discarded [DEFAULT: "+str(MINIMUM_DEPTH)+"].", required=False, default=MAXIMUM_GNOMAD)
    parser.add_argument('-gerp','--gerp',action='store',type=float,help="Minimum GERP++ score for the variant to be considered deleteriuos [DEFAULT: "+str(GERP_minimum)+"].", required=False, default=GERP_minimum)
    parser.add_argument('-s','--score',action='store',type=float,help="Minimum number of pathogenic score predictors to be evaluated as 'DELETERIOUS' to accept a variant with an uncertain ClinVar score [DEFAULT: "+str(MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS)+"].", required=False, default=MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS)
    
    ### OUTPUT - File
    parser.add_argument('-o','--outputpass',action='store',type=str,help="Path of the VCF output data filtered by PASS.", required=False, default="output_passing.vcf")
    arg=parser.parse_args()

    INPUT=arg.input
    
    ### GENETIC Regions form BED files and UCSC Genome Browser
    GENE_SEARCHED = getGeneticRegions(arg.genelist, arg.bedgenefile)
    
    ### FILTER IMPACT - Parameters
    FILTER_IMPACT=arg.filterimpact
    if FILTER_IMPACT is not None :
        FILTER_IMPACT=set(FILTER_IMPACT.upper().split(','))
        print("IMPACT FILTER: {:}".format(", ".join(FILTER_IMPACT)))
    
    ############################################################### FILTER
    def doAF_bigger(value):
        if isinstance(value,tuple) :
            value=value[-1]
        return value is None or value == "." or float(value) >= arg.minAF
    
    def doGNOMAD_rare(value):
        return value is None or value == "." or float(value) <= arg.gnomad
    
    def doDEPTH(value):
        return value is not None and value>arg.minDP
    
    dictionaryFilter={}
    dictionaryFilter['FORMAT.AF']=doAF_bigger
    dictionaryFilter[arg.population]=doGNOMAD_rare
    dictionaryFilter['DP']=doDEPTH
    
    print("INPUT FILE: {:}".format(INPUT))
    reader=readerVCF_AnnotationFilter_Report(INPUT, selectedGenes=GENE_SEARCHED, outPass=arg.outputpass, filteringCondition=dictionaryFilter, filteringPriority=None)
    reader.selectedImpactSet=FILTER_IMPACT
    reader.MINIMUM_REQUIRED_PATHOGENIC_PREDICTORS=arg.score
    
    ### SET the parameters - ERROR in the VCF READING are BLOCKING
    reader.BLOCKING_ERROR=arg.set_blocking_errors
    ### SET the parameters - INCLUDE CONTIGs in the HEADER
    reader.INCLUDE_contig=arg.set_include_contigs
    
    ### READ the VCF Data
    start_time = time.time()
    counterRow=reader.loadData()
    print("VCF [{:d}] - Loading time {:.3f} s".format(counterRow,time.time() - start_time))
    reader.close()
    
    ### PRINT a REPORT of the VCF
    print("Type: %s"%reader.vcfType)
    print("# Dictionary FORMAT")
    print(reader.dictionaryFormat_toString())
    print()
    print("# Dictionary INFO")
    print(reader.dictionaryInfo_toString())
    print()
    
    if reader.selectedImpactSet is not None :
        print("Filter applied to variant Impact: {:}".format(" ".join(reader.selectedImpactSet)))
    
    print(reader.headerColumns)
    if counterRow>0 :
        print("### Passing VARIANTS: {:d}".format(reader.counterPassing))   
    else:
        print("### NON PASSING VARIANTS - Not Found")    

    print()
    