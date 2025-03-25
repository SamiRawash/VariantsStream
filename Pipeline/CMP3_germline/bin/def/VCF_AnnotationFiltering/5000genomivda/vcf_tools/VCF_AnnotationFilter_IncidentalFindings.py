#! /usr/bin/env python3
'''
Created on Jan 16, 2023

@author: flanduzzi
'''  
import argparse, time
from vcf_tools.VCF_AnnotationFilter import readerVCF_AnnotationFilter_BaseClass, getGeneticRegions
DEFAULT_BED="/opt/BED/ACMG_IncidentalFindings.bed" # None  # 

CLINVAR_pathogenic="pathogenic"
CLINVAR_conflicting="conflicting_interpretations_of_pathogenicity"

class readerVCF_AnnotationFilter_Report(readerVCF_AnnotationFilter_BaseClass):
    
    def __init__(self, source, selectedGenes=None, outPass="output_passing.vcf"):
        super().__init__(source, selectedGenes=selectedGenes, outPass=outPass)
    
    def doWithPassingAnnotatedVariant(self,variant):
        clinvar=variant.info['CLNSIG']
        if CLINVAR_pathogenic in clinvar.lower() and CLINVAR_conflicting not in clinvar.lower() :
            self.outfilePass.write(variant.txt + '\n')
            self.counterPassing += 1
            return True
        return False
    
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
    parser.add_argument('-bed','--bedgenefile',action='store',type=str,help=f"BED file containing the list of genes to be searched [DEFAULT: ACMG v3.1 - location {DEFAULT_BED}].", required=False, default=DEFAULT_BED)
    
    ### OUTPUT - File
    parser.add_argument('-o','--outputpass',action='store',type=str,help="Path of the VCF output data filtered by PASS.", required=False, default="output_passing.vcf")
    arg=parser.parse_args()

    INPUT=arg.input
    
    ### GENETIC Regions form BED files and UCSC Genome Browser
    GENE_SEARCHED = getGeneticRegions(arg.genelist, arg.bedgenefile)
    
    print("INPUT FILE: {:}".format(INPUT))
    reader=readerVCF_AnnotationFilter_Report(INPUT, selectedGenes=GENE_SEARCHED, outPass=arg.outputpass)
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
    
    print(reader.headerColumns)
    if counterRow>0 :
        print("### Passing VARIANTS: {:d}".format(reader.counterPassing))   
    else:
        print("### Passing VARIANTS - Not Found")    

    print()
    