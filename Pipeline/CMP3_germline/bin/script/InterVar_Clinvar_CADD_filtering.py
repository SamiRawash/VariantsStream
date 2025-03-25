#! /usr/local/bin/python3
import argparse, time

def split_line_by_separator(vcf_line,separator="\t"):
    splitted_row=vcf_line.split(separator)
    return splitted_row


def get_InterVar_from_vcf_line(element):
    if element=="Benign" or element=="Likely_benign":
        return 0
    elif element=="Pathogenic" or element=="Likely_pathogenic":
        return 1
    elif element==".":
        return 2
    else :
        return None

def get_Clinvar_from_vcf_line(element):
    return element
   
def get_ClinStars_from_vcf_line(element):
    if "no_assertion" in element:
        return 0
    elif "single_submitter" in element or "conflicting" in element:
        return 1
    elif "multiple" in element:
        return 2
    elif "reviewed" in element:
        return 3
    elif "practice_guideline" in element:
        return 4
    else:         
        return -1

def get_AF_from_vcf_line(element):
    if element == "." or element == "" :
        return None
    else :
        return float(element)
                   
def get_CADD_from_vcf_line(element):
    if element == ".":
        return None
    else :
        return float(element)                           

dictionary_filters={
    "InterVar_automated" : get_InterVar_from_vcf_line ,
    "CLNSIG" : get_Clinvar_from_vcf_line , #lambda x : x.strip() if x is not None else None ,
    "CLNREVSTAT" : get_ClinStars_from_vcf_line ,
    "AF" : get_AF_from_vcf_line ,
    "CADD_phred" : get_CADD_from_vcf_line
    }
'''
def update_results(InterVar_val,Clinvar,Stars,CADD,FREQ_AF, results):
    if (InterVar_val,Clinvar,Stars,CADD,FREQ_AF) in results :
        results[(InterVar_val,Clinvar,Stars,CADD,FREQ_AF)]+=1
    else :
        results[(InterVar_val,Clinvar,Stars,CADD,FREQ_AF)]=1
    return results
'''
if __name__=="__main__":
    parser=argparse.ArgumentParser(description="Get gene from VCF")
    parser.add_argument('-v','--vcf',action='store',type=str,help='Path to the VCF file',required=True,default=None)
    parser.add_argument('-o','--output',action='store',type=str,help='Path to the output VCF file',required=False,default="Filtered.vcf")
    parser.add_argument('-wm','--writing_mode',action='store',type=str,help='Writing mode accepted values are: write ("w"), append ("a") [DEFAULT: w] ',required=False,default="w")
    args=parser.parse_args()
    vcf_file_path=args.vcf
    output_file=args.output
    clinvar_stats_dict={'no_assertion':'0','single_submitter':'1','conflicting':'1','multiple':'2','reviewed':'3','practice_guideline':'4'}
    
    dictionary_results={}
    start_time = time.time()
    with open(output_file, args.writing_mode) as vcf_output :
        with open(vcf_file_path, "r") as vcf_file :
            line=vcf_file.readline()
            n=0  
            m=0
            while line!="":
                if line[0]=="#":
                    vcf_output.write(line)             
                else:
                    m+=1
                    splitted_info=line.strip().split("\t")[7].split(";")
                    #info={}
                    InterVar_val=None
                    Clinvar=None
                    Stars=None
                    CADD=None
                    FREQ_AF=None
                    for element in splitted_info :
                        el_vector=element.split("=", maxsplit=1)
                        if len(el_vector)>0 and el_vector[0] in dictionary_filters.keys() :
                            #info[el_vector[0]]=dictionary_filters[el_vector[0]](el_vector[1])
                            if el_vector[0] == "InterVar_automated" :
                                InterVar_val=get_InterVar_from_vcf_line(el_vector[1])
                            elif el_vector[0] == "CLNSIG" :
                                Clinvar=get_Clinvar_from_vcf_line(el_vector[1])
                            elif el_vector[0] == "CLNREVSTAT" :
                                Stars=get_ClinStars_from_vcf_line(el_vector[1])
                            elif el_vector[0] == "CADD_phred" :
                                CADD=get_CADD_from_vcf_line(el_vector[1])
                            elif el_vector[0] == "AF" :
                                FREQ_AF=get_AF_from_vcf_line(el_vector[1])
                    
                    if InterVar_val is not None and InterVar_val==1:
                        vcf_output.write(line)
                        #dictionary_results=update_results(InterVar_val,Clinvar,Stars,CADD,FREQ_AF,dictionary_results)
                        n+=1
                    elif InterVar_val is not None and InterVar_val==0 :
                        pass
                    else : # A.K.A. InterVal_val=2
                        if Clinvar is not None and (Clinvar=="Benign" or Clinvar=="Likely_benign" or Clinvar=="Benign/Likely_benign") and Stars>=1:    
                            pass
                        elif Clinvar is not None and (Clinvar=="Pathogenic" or Clinvar=="Likely_pathogenic" or Clinvar=="Pathogenic/Likely_pathogenic") and Stars>=1:    
                                vcf_output.write(line)
                                #dictionary_results=update_results(InterVar_val,Clinvar,Stars,CADD,FREQ_AF,dictionary_results)
                                n+=1
                        else :
                            if CADD is not None and CADD >=20 :    
                                vcf_output.write(line)
                                #dictionary_results=update_results(InterVar_val,Clinvar,Stars,CADD,FREQ_AF,dictionary_results)
                                n+=1
                            elif CADD is None and (FREQ_AF is None or FREQ_AF<0.01):
                                vcf_output.write(line)
                                #dictionary_results=update_results(InterVar_val,Clinvar,Stars,CADD,FREQ_AF,dictionary_results)
                                n+=1
                                                                                                      
                line=vcf_file.readline()
        
    # The line below is the original one, commented out to add exception if division by zero
    #print(f"Total Variants [{m}] Surviving Variants[{n}] : {(float(n)/m*100):.5}% - Execution time {time.time()-start_time:.3f} s")    
    if m != 0:
        print(f"Total Variants [{m}] Surviving Variants[{n}] : {(float(n)/m*100):.5}% - Execution time {time.time()-start_time:.3f} s")
    else:
        print(f"Total Variants [{m}] Surviving Variants[{n}] ((float(n)/m*100):.5% not computed to avoid division by zero) - Execution time {time.time()-start_time:.3f} s")
    
    
    
    #keys=list(dictionary_results.keys())
    #keys.sort()
    #for InterVar_val,Clinvar,Stars,CADD,FREQ_AF in keys :
    #    val = dictionary_results[(InterVar_val,Clinvar,Stars,CADD,FREQ_AF)]
    #    print(f"{val}\t{InterVar_val}\t{Clinvar}\t{Stars}\t{CADD}\t{FREQ_AF}")

