#! /usr/bin/env python3
'''
Created on Jun 18, 2021

@author: flanduzzi
'''  
import argparse, time
from vcf_tools.VCF_AnnotationFilter import readerVCF_AnnotationFilter, getGeneticRegions

VARIANT_IMPACT="HIGH,MODERATE"
REFERENCE_POPULATION='gnomAD_genome_ALL'
MAXIMUM_ALLEL_FREQUENCY=0.05

################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="FILTER ANNOTATED VCF: considering a VCF annotated  with PASS, SNPEFF, ANNOVAR the script select the passing variants that refer to a preselected GENES_SET. The set of variants is the filter again by the IMPACT and AllelFrequency in a determined population. ")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file (produced by SPNEff, ANNOVAR or VEP).", required=True, default=None)
    parser.add_argument('-blocking','--set_blocking_errors',action='store_true', help="Flag forcing the end of the execution with the raise of an exception in case of an error during the VCF readings [DEFAULT: False].", default=False)
    parser.add_argument('-contigs','--set_include_contigs',action='store_true', help="Flag allowing the inclusion of the CONTIG raws in the VCF header [DEFAULT: False].", default=False)
    
    ### PARAMETERs
    parser.add_argument('-g','--genelist',action='store',type=str,help="List of genes to be searched separeted by commas (no spaces).", required=False, default=None)
    parser.add_argument('-bed','--bedgenefile',action='store',type=str,help="BED file containing the list of genes to be searched.", required=False, default=None)
    parser.add_argument('-f','--filterimpact',action='store',type=str,help="String containing the variant gene impact separated by commas [DEFAULT: "+VARIANT_IMPACT+"].", required=False, default=VARIANT_IMPACT)
    parser.add_argument('-p','--population',action='store',type=str,help="String describing the field from which extract REFERENCE POPULATION for the allel frequency (AF). Notice that if the information is stored in the 'CSQ' field (a.k.a. VEP annotation) only the population ALL is checked. [DEFAULT: "+REFERENCE_POPULATION+"].", required=False, default=REFERENCE_POPULATION)
    parser.add_argument('-af','--maxAF',action='store',type=float,help="Set the maximum ALLEL FREQUENCY in the reference POPULATION above which discard the variant [DEFAULT: "+str(MAXIMUM_ALLEL_FREQUENCY)+"].", required=False, default=MAXIMUM_ALLEL_FREQUENCY)
    
    ### OUTPUT - File        
    parser.add_argument('-op','--outputpass',action='store',type=str,help="Path of the VCF output data filtered by PASS.", required=False, default="output_passing.vcf")
    parser.add_argument('-oi','--outputimpact',action='store',type=str,help="Path of the VCF output data filtered by PASS and IMPACT.", required=False, default="output_passing_HighModerate.vcf")
    parser.add_argument('-oa','--outputaf',action='store',type=str,help="Path of the VCF output data filtered by PASS, IMPACT, COSMIC and Allel Frequency.", required=False, default="output_passing_HighModerate_Cosmic_AF.vcf")
    arg=parser.parse_args()

    INPUT=arg.input
    
    ### GENETIC Regions form BED files and UCSC Genome Browser
    GENE_SEARCHED = getGeneticRegions(arg.genelist, arg.bedgenefile)
    
    ### FILTER IMPACT - Parameters
    FILTER_IMPACT=arg.filterimpact
    if FILTER_IMPACT is not None :
        FILTER_IMPACT=set(FILTER_IMPACT.upper().split(','))
        print("IMPACT FILTER: {:}".format(", ".join(FILTER_IMPACT)))
    
    
    print("INPUT FILE: {:}".format(INPUT))
    reader=readerVCF_AnnotationFilter(INPUT, selectedGenes=GENE_SEARCHED, setImpact=FILTER_IMPACT, outPass=arg.outputpass, outImpact=arg.outputimpact, outAF=arg.outputaf)
    ### SET the parameters for the ALLEL FREQUENCY in the POPULATION
    
    ### SET the parameters for the ALLEL FREQUENCY in the POPULATION
    REF_POPULATION=arg.population
    reader.maxAF=arg.maxAF
    reader.field_AF=REF_POPULATION
    
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
        print("### Impact VARIANTS: {:d}".format(reader.counterHighModerate))   
        print("### AllelFrequency VARIANTS reference population {:}<{:.3f}: {:d}".format(reader.field_AF, reader.maxAF, reader.counterAllelFequency))   

    else:
        print("### NON PASSING VARIANTS - Not Found")    

    print()
    