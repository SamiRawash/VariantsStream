#! /usr/bin/env python3
import argparse, time, copy
from vcf_tools.VCF_AnnotationFilter import readerVCF_AnnotationFilter_BaseClass
from vcf_tools.cancervar.CancerVar_Criteria import initialize_all_database,assign,CBP_TRESHOLD_LIKELY_PATHOGENIC_MIN

BASE_FOLDER=""
DEFAULT_CANCER_TYPE="All Tumors"
DEFAULT_DATABASE_LOF=BASE_FOLDER+"data/LOF.genes.exac_me_cancers"
DEFAULT_DATABASE_MARKER=BASE_FOLDER+"data/cancervar.out.txt"
DEFAULT_DATABASE_PATHWAYS=BASE_FOLDER+"data/cancers_genes.list_kegg.txt"
DEFAULT_DATABASE_CANCER_GENES=BASE_FOLDER+"data/cancer_census.genes"
DEFAULT_DATABASE_CANCER_TYPES=BASE_FOLDER+"data/cancervar.cancer.types"

class VCF_addCancerVar(readerVCF_AnnotationFilter_BaseClass):
    
    pathogenic_variants=[]
    
    def doWithPassingAnnotatedVariant(self, variant):
        
        tier,score,cbp,Therap_list,Diag_list,Prog_list=assign(variant, self.cancer_type)
        #tier,score,cbp=assign(variant, self.cancer_type)
        if score>= CBP_TRESHOLD_LIKELY_PATHOGENIC_MIN: 
            v=copy.deepcopy(variant)
            self.pathogenic_variants.append( (v,tier,score,cbp,Therap_list,Diag_list,Prog_list) )
            

        new_row=variant.txt.split("\t")
        new_row[7]+=f";CancerVar={tier};CancerVar_score={score};EVS={cbp}"
        self.outfilePass.write("\t".join(new_row) + '\n')
        self.counterPassing += 1
    
    def doWithHeaderColumn(self, line):
        self.headerTxt=self.headerTxt[:-(len(line)+1)]
        self.headerTxt+='##INFO=<ID=CancerVar,Number=.,Type=String,Description="CancerVar Tier accordingly to the Assoc.MedicalPathology/AmericanSocietyClinicalOncology/CollegeAmericanPathologists 2017 guidelines">\n'
        self.headerTxt+='##INFO=<ID=CancerVar_score,Number=.,Type=Integer,Description="CancerVar Score accumulated in after evaluation of the 12 clinical-based evidence.">\n'
        self.headerTxt+='##INFO=<ID=EVS,Number=.,Type=String,Description="Score in the 12 clinical-based evidence to predict the clinical significance of somatic variants.">\n'
        self.headerTxt+=line+"\n"
        return readerVCF_AnnotationFilter_BaseClass.doWithHeaderColumn(self, line)
        
################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(prog="vcf_cancervar", description="Using a VCF annotated with ANNOVAR (vers.) evaluate the 12 clinical-based evidence to predict the clinical significance of somatic variants accordingly to the Assoc.MedicalPathology/AmericanSocietyClinicalOncology/CollegeAmericanPathologists 2017 guidelines. ")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file annotated with ANNOVAR.", required=True, default=None)
    parser.add_argument('-c','--cancer_type',action='store',type=str,help=f"Type of cancer under analysis [DEFAULT:{DEFAULT_CANCER_TYPE}].", required=False, default=DEFAULT_CANCER_TYPE)
    parser.add_argument('-blocking','--set_blocking_errors',action='store_true', help="Flag forcing the end of the execution with the raise of an exception in case of an error during the VCF readings [DEFAULT: False].", default=False)
    
    parser.add_argument('-lof','--lof',action='store',type=str,help=f"Database of the LoF genes [DEFAULT: {DEFAULT_DATABASE_LOF}].", required=False, default=DEFAULT_DATABASE_LOF)
    parser.add_argument('-mk','--marker',action='store',type=str,help=f"Database of the marker genes [DEFAULT: {DEFAULT_DATABASE_MARKER}].", required=False, default=DEFAULT_DATABASE_MARKER)
    parser.add_argument('-pt','--pathways',action='store',type=str,help=f"Database of the pathways [DEFAULT: {DEFAULT_DATABASE_PATHWAYS}].", required=False, default=DEFAULT_DATABASE_PATHWAYS)
    parser.add_argument('-cgc','--cancergenes',action='store',type=str,help=f"Database of the cancer genes [DEFAULT: {DEFAULT_DATABASE_CANCER_GENES}].", required=False, default=DEFAULT_DATABASE_CANCER_GENES)
    parser.add_argument('-ct','--cancertypes',action='store',type=str,help=f"Database of the cancer types [DEFAULT: {DEFAULT_DATABASE_CANCER_TYPES}].", required=False, default=DEFAULT_DATABASE_CANCER_TYPES)
    parser.add_argument('-ev','--userevidence',action='store',type=str,help=f"Database of the user defined evidence [DEFAULT: None].", required=False, default=None)
    
    ### PARAMETERs
    parser.add_argument('-o','--output',action='store',type=str,help="Output file path used to store the tabular results .", required=False, default="Table_output.vcf")
    arg=parser.parse_args()

    INPUT=arg.input
    #OUTPUT=arg.output
    output=arg.output
    
    initialize_all_database(arg.marker, arg.lof, arg.pathways, arg.cancergenes, arg.cancertypes)
    #cancervar_d,cancervar_markers_dict=loadCancervarDB(arg.marker)
    #lof_genes_dict=loadLofDB(arg.lof)
    #cancer_pathway_dict=loadCancerPathways(arg.pathways)
    #cancers_gene_dict=loadCancerGenes(arg.cancergenes)
    #cancer_types_dict=loadCancerTypes(arg.cancertypes)
    
    print("INPUT FILE: {:}".format(INPUT))
    reader=VCF_addCancerVar(INPUT, outPass=output)
    reader.cancer_type=arg.cancer_type
    ### SET the parameters - ERROR in the VCF READING are BLOCKING
    reader.BLOCKING_ERROR=arg.set_blocking_errors### Loading of the VCF Data
    
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
    if len(reader.pathogenic_variants)>0 :
        print(f"Cancervar Pathogenic/LikelyPathogenic variants [{len(reader.pathogenic_variants)}/{reader.counterPassing}]:")
        for variant,tier,score,cbp,Therap_list,Diag_list,Prog_list in reader.pathogenic_variants :
            print(f"{variant.info['Gene.refGene_latest']:15}\t{variant.chromosome}:{variant.position}:{variant.reference}:{variant.alteration:<6}\t{tier:<20}\tSCORE:{score}\tEVS={cbp}\tTherapy:{Therap_list}\tDiagnosis:{Diag_list}\tPrognosis:{Prog_list}")
    else:
        print(f"Cancervar Pathogenic/LikelyPathogenic variants [{len(reader.pathogenic_variants)}/{reader.counterPassing}]: None")
    
    
    
        
            