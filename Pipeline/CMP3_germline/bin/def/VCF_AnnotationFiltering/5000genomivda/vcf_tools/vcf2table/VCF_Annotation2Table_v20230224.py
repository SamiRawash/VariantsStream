#! /usr/bin/env python3
'''
Created on Feb 24, 2023

@author: flanduzzi
'''  
import argparse, time
from vcf_tools.parserAnnotationVCF import parserAnnotator_ANN,parserAnnotator_CSQ,parserAnnotator_LOF
from vcf_tools.TypePolymorfism import typePolymorfism_Extended
from vcf_tools.CommonTools import chiaveVal, chromosomeCode, chromosomeCode2sting, templateCompressReader, DICTIONARY_GENOTYPE
from vcf_tools.TypeAnnotation import typeInfo
from vcf_tools.vcf2table.varsome_criteria import loadVarsomeCriteriaTable

BASE_FOLDER=""
DEFAULT_VARSOME_CRITERIA_TABLE=BASE_FOLDER+"data/Varsome_CriteriaOnPathSorePredictors.txt"
DEFAULT_SNPEFF_INTERESTING_IMPACT=set(["HIGH","MODERATE"])
DEFAULT_ANNOVAR_GNOMAD_FIELDS=set(["AF","AF_raw","AF_male","AF_female","AF_afr","AF_ami","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth","AF_sas"])
EMPTY_FIELD="."
OMIM_DOMINANCE_DICT={"Autosomal dominant" : "AD", "Autosomal recessive": "AR"}

#PATHOGENIC_SCORE_PROTEIN_FUNCTION=set(["SIFT_pred", "SIFT4G_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred", "LRT_pred", "MutationAssessor_pred", "FATHMM_pred", "PROVEAN_pred", "MetaSVM_pred", "MetaLR_pred", "MetaRNN_pred"])
#PATHOGENIC_SCORE_CLINICAL=set(["MutationTaster_pred",])

PATHOGENIC_SCORE_LIST_INDIVIDUAL=['CADD_phred','DANN_score','DEOGEN2_score','Eigen-raw_coding','Eigen-PC-raw_coding',
                  'FATHMM_score','fathmm-MKL_coding_score','fathmm-XF_coding_score','LIST-S2_score','LRT_score',
                  'MutPred_score','MutationAssessor_score','MutationTaster_score','PROVEAN_score','Polyphen2_HDIV_score','Polyphen2_HVAR_score','PrimateAI_score',
                  'SIFT_score','SIFT4G_score','dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE','MVP_score']
PATHOGENIC_SCORE_LIST_METASCORE=['BayesDel_addAF_score','BayesDel_noAF_score','MetaLR_score','MetaRNN_score','MetaSVM_score','REVEL_score','M-CAP_score']
PATHOGENIC_SCORE_LIST_CONSERVATION=['phastCons100way_vertebrate','phyloP30way_vertebrate','SiPhy_29way_logOdds','GERP++_NR','GERP++_RS']

def convertLoF2string(field) :
    return "{:}".format(field)

def convertNmd2string(field) :
    return "{:}".format(field)

def convertClnrevstat2star(field):
    if "no_assertion" in field:
        return 0
    elif "single_submitter" in field or "conflicting" in field:
        return 1
    elif "multiple" in field:
        return 2
    elif "reviewed" in field:
        return 3
    elif "practice_guideline" in field:
        return 4
    else:         
        return -1

def getGnomAD_MaxPopulationFrequency(fields_gnomAD):
#    values=set([float(x) if x is not None and x != "." else None for x in fields_gnomAD])
#    return max( values ) if len(values)>0 and next(iter(values)) is not None else EMPTY_FIELD
    # Convert values to float, replace None for "." or missing values
    values = {float(x) for x in fields_gnomAD if x is not None and x != "."}
    
    # Check if the set is not empty before trying to get the max value
    return max(values) if values else EMPTY_FIELD




def getOmimFields(omim_field):
    ### Ex. bed2=Name\x3dInsensitivity_to_pain,_congenital,_with_anhidrosis,_256800_(3),_Autosomal_recessive
    if omim_field is not None and omim_field.startswith("Name") :
        omim=""
        omim_inheritance=""
        omim_pathology=""
        #val= omim_field.decode('utf-8')#decode("utf-8", "replace")
        #val=val.decode('utf-8')
        #str(omim_field, encoding='utf-16') #.decode('utf-16le', errors='ignore').encode('utf8')
        modified=omim_field[8:].replace('\\x',";")
        modified=modified.replace(";3b",";3d")
        modified=modified.replace(";3d",";")
        try: 
            for row in modified.split(';'):
                splitted=row.replace("_"," ").split(",")
                if len(splitted)>2 :
                    val_omim=["","",""]
                    position=0
                    count=0
                    while position<3 and count<len(splitted):
                        tmp_omim=splitted[count].strip()
                        count+=1
                        if ("{" in tmp_omim and ("}" not in tmp_omim) )  or ("]" in tmp_omim and ("[" not in tmp_omim)) :
                            #tmp_omim=tmp_omim.replace("{","")
                            while (("}" not in splitted[count]) and ("]" not in splitted[count])) and count<len(splitted):
                                tmp_omim+="/"+ splitted[count].strip()
                                count+=1
                            tmp_omim+="/"+ splitted[count].strip()
                        val_omim[position]+=tmp_omim
                        position+=1
                    omim_inh_tmp=val_omim[2].strip()
                    omim+=val_omim[1].strip() +"," #f"https://www.omim.org/entry/{val_omim[1].strip()}"+","
                    omim_inheritance+=(OMIM_DOMINANCE_DICT[omim_inh_tmp] if omim_inh_tmp in OMIM_DOMINANCE_DICT else omim_inh_tmp) +","
                    omim_pathology+=val_omim[0].strip()+","
                #print(f"{splitted}\t\t{omim};{omim_inheritance};{omim_pathology}")
        except Exception :
            print(f"Error: {omim_field}")
        return omim[:-1],omim_inheritance[:-1],omim_pathology[:-1]
        
    else :
        return (EMPTY_FIELD,EMPTY_FIELD,EMPTY_FIELD)

def getOmimGeneCode(omim_field):
    ### Ex. bed3=Name\x3d191315
    if omim_field is not None and omim_field.startswith("Name") :
        rpl=omim_field[8:]
        try:
            return int(rpl)
        except:
            return rpl
    return EMPTY_FIELD

def getFORMAT_2String(formatColumn, index=None):
    DP=EMPTY_FIELD
    AD_ref=EMPTY_FIELD
    AD_alt=EMPTY_FIELD
    GT=EMPTY_FIELD
    GQ=EMPTY_FIELD
    
    if index is not None :
        try :
            DP=int(formatColumn['DP'][index])
        except :
            pass
        
        if 'AD' in formatColumn :
            val=formatColumn['AD'][index].split(",")
            AD_ref=EMPTY_FIELD
            try:
                AD_ref=int(val[0]) if val[0]!="." else EMPTY_FIELD
            except: pass
            
            AD_alt=EMPTY_FIELD
            try:
                AD_alt=int(val[1]) if len(val)>1 and val[1]!="." else EMPTY_FIELD 
            except: pass
        
        GT=EMPTY_FIELD if 'GT' not in formatColumn else formatColumn['GT'][index]
        if GT in DICTIONARY_GENOTYPE :
            GT=DICTIONARY_GENOTYPE[GT]
    
        GQ=EMPTY_FIELD if 'GQ' not in formatColumn else formatColumn['GQ'][index]
        if GQ in DICTIONARY_GENOTYPE :
            GQ=DICTIONARY_GENOTYPE[GQ]
        
    else :
        try :
            DP=int(formatColumn['DP'])
        except :
            pass
        
        if 'AD' in formatColumn :
            val=formatColumn['AD'].split(",")
            AD_ref=int(val[0]) if val[0]!="." else EMPTY_FIELD
            AD_alt=int(val[1]) if len(val)>1 and val[1]!="." else EMPTY_FIELD
        
        GT=EMPTY_FIELD if 'GT' not in formatColumn else formatColumn['GT']
        if GT in DICTIONARY_GENOTYPE :
            GT=DICTIONARY_GENOTYPE[GT]
    
        GQ=EMPTY_FIELD if 'GQ' not in formatColumn else formatColumn['GQ']
        if GQ in DICTIONARY_GENOTYPE :
            GQ=DICTIONARY_GENOTYPE[GQ]
        
    return DP,AD_ref,AD_alt,GT,GQ
        
def loadCSV_Column(source_path, separator="\t", column_key=0, column_value=1, mode="r", header=True):
    set_of_names={}
    with open(source_path, mode=mode) as csv:
        if header :
            csv.readline()
        for line in csv.readlines() :
            if line[0]!="#" and line.strip()!="" :
                splitted=line.strip().split(separator)
                key=splitted[column_key].strip()
                set_of_names[key]=splitted[column_value].strip()
    return set_of_names


class readerVCF_AnnotationFilter(templateCompressReader):
    
    varsomeCriteria=None
    
    def __init__(self, source, selectedfieldsINFO=None, output="Table.csv"):
        super().__init__()
        self.source_path=source
        self.vcfType="unknown"
        self.headerColumns=""
            
        ### Internal Variables
        self.TEMPLATE_INFO=""
        self.outfile=open(output, 'w')
        self.DUMP_SAMPLE_NAME=False
        
        self.KEEP_ONLY_PASS=False
        self.KEEP_ONLY_TRANSCRIPT_WITH_IMPACT=False
        self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT=False
        
        self.SNPEFF_INTERESTING_IMPACT=DEFAULT_SNPEFF_INTERESTING_IMPACT
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
        
        ### LIST of GENE NAMEs
        self._list_of_genes=None
        self._list_of_canonical_transcript=None
    
    ###### Function for LOADING External INFORMATION FILEs ################################################################
    def setGeneNames(self, gene_names):
        if isinstance(gene_names, set) :
            self._list_of_genes=gene_names
        else :
            raise Exception("GeneNames must be a SET")
    
    def setCanonicalTranscript(self, can_trans):
        if isinstance(can_trans, set) :
            self._list_of_canonical_transcript=can_trans
        else :
            raise Exception("CanonicalTranscript must be a SET")
    ########################################################################################################################

    def cleanListOfInterestingGenes(self):
        self._list_of_genes.clear()
    
    def doWithLines(self, line, variantsCounter):
        ### Increase the variants counter
        variantsCounter = variantsCounter + 1
        variant, filtro = self.line2typePolymorfism(line)
        
        if self.KEEP_ONLY_PASS and filtro!='PASS':
            return variantsCounter
        
        if self.KEEP_ONLY_TRANSCRIPT_WITH_IMPACT :
            self.keepOnlyImpactTranscripts_SnpEff(variant)
            
        self.markSnpeffTranscript(variant, keepOnlyCanonical=self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT)
        
        info=variant.info
        if 'ANN' in info and len(info['ANN'])>0:
            
            flag_lof=True if 'LOF' in info else False
            lof=convertLoF2string(info['LOF']) if flag_lof else EMPTY_FIELD
            flag_nmd=True if 'NMD' in info else False
            nmd=convertNmd2string(info['NMD']) if flag_nmd else EMPTY_FIELD
            
            chromosome=chromosomeCode2sting(variant.chromosome)
            position=variant.position
            ref=variant.reference
            alt=variant.alteration
            
            avsnp=info['avsnp_latest'] if 'avsnp_latest' in info else EMPTY_FIELD
            annovar_FuncRefGene=info['Func.refGene_latest'] if 'Func.refGene_latest' in info else EMPTY_FIELD
            annovar_ExFunc=info['ExonicFunc.refGene_latest'] if 'ExonicFunc.refGene_latest' in info else EMPTY_FIELD
            gnomAD=getGnomAD_MaxPopulationFrequency([info[gnmkey] if gnmkey in info else None for gnmkey in DEFAULT_ANNOVAR_GNOMAD_FIELDS])
            omim,omim_inheritance,omim_pathology=getOmimFields(info['bed2']) if 'bed2' in info else (EMPTY_FIELD,EMPTY_FIELD,EMPTY_FIELD)
            omim_genecode=getOmimGeneCode(info['bed3']) if 'bed3' in info else EMPTY_FIELD
            ACMG=EMPTY_FIELD
            if  self.columnFormat_SampleNumber == 2:
                ACMG=info['CancerVar'] if 'CancerVar' in info else EMPTY_FIELD
                ACMG_score=info['CancerVar_score'] if 'CancerVar_score' in info else EMPTY_FIELD
                ACMG=f"{ACMG} [{ACMG_score}]"
            else:
                ACMG=info['InterVar_automated'] if 'InterVar_automated' in info else EMPTY_FIELD 
            Clinvar=info['CLNSIG'] if 'CLNSIG' in info else EMPTY_FIELD
            Clinvar_star=convertClnrevstat2star(info['CLNREVSTAT']) if 'CLNREVSTAT' in info else EMPTY_FIELD
            Clinvar_info=info['CLNDN'] if 'CLNDN' in info else EMPTY_FIELD
            cadd=info['CADD_phred'] if 'CADD_phred' in info else EMPTY_FIELD
            annovar_NM=info['AAChange.refGene_latest'] if 'AAChange.refGene_latest' in info else EMPTY_FIELD
            cytoBand=info['cytoBand'] if 'cytoBand' in info else EMPTY_FIELD
            
            ##############################
            ### INDIVIDUAL SCORE PREDICTOR
            ##############################
            score_pred_individual_scor,score_pred_individual_pat,score_pred_individual_ben,INDIVIDUAL_SCOREPREDICTOR_BODY=self.aggregateScorePredictor(info, PATHOGENIC_SCORE_LIST_INDIVIDUAL)
            """
            cadd=info['CADD_phred'] if 'CADD_phred' in info else EMPTY_FIELD
            cadd_pred=self.varsomeCriteria['CADD_phred'][4](cadd)
            META_SCOREPREDICTOR_HEAD=""
            META_SCOREPREDICTOR_BODY=""

            dann=info['DANN_score'] if 'DANN_score' in info else EMPTY_FIELD
            dann_pred=self.varsomeCriteria['DANN_score'][4](cadd)
            
            deogen2=info['DEOGEN2_score'] if 'DEOGEN2_score' in info else EMPTY_FIELD
            deogen2_pred=self.varsomeCriteria['DEOGEN2_score'][4](deogen2)
            
            eigen=info['Eigen-raw_coding'] if 'Eigen-raw_coding' in info else EMPTY_FIELD
            eigen_pred=self.varsomeCriteria['Eigen-raw_coding'][4](eigen)
            
            eigenpc=info['Eigen-PC-raw_coding'] if 'Eigen-PC-raw_coding' in info else EMPTY_FIELD
            eigenpc_pred=self.varsomeCriteria['Eigen-PC-raw_coding'][4](eigenpc)
            
            fathmm=info['FATHMM_score'] if 'FATHMM_score' in info else EMPTY_FIELD
            fathmm_pred=self.varsomeCriteria['FATHMM_score'][4](fathmm)
            
            fathmmmkl=info['fathmm-MKL_coding_score'] if 'fathmm-MKL_coding_score' in info else EMPTY_FIELD
            fathmmmkl_pred=self.varsomeCriteria['fathmm-MKL_coding_score'][4](fathmmmkl)
            
            fathmmxf=info['fathmm-XF_coding_score'] if 'fathmm-XF_coding_score' in info else EMPTY_FIELD
            fathmmxf_pred=self.varsomeCriteria['fathmm-XF_coding_score'][4](fathmmxf)
            
            lists2=info['LIST-S2_score'] if 'LIST-S2_score' in info else EMPTY_FIELD
            lists2_pred=self.varsomeCriteria['LIST-S2_score'][4](lists2)
            
            lrt=info['LRT_score'] if 'LRT_score' in info else EMPTY_FIELD
            lrt_pred=self.varsomeCriteria['LRT_score'][4](lrt)
            
            mutpred=info['MutPred_score'] if 'MutPred_score' in info else EMPTY_FIELD
            mutpred_pred=self.varsomeCriteria['MutPred_score'][4](mutpred)
            mutass=info['MutationAssessor_score'] if 'MutationAssessor_score' in info else EMPTY_FIELD
            mutass_pred=self.varsomeCriteria['MutationAssessor_score'][4](mutass)
            muttas=info['MutationTaster_score'] if 'MutationTaster_score' in info else EMPTY_FIELD
            muttas_pred=self.varsomeCriteria['MutationTaster_score'][4](muttas)
            #META_SCOREPREDICTOR_HEAD+=f"Polyphen2_HDIV\tPolyphen2_HDIV_pred\Polyphen2_HVAR\tPolyphen2_HVAR_pred\tPrimateAI\tPrimateAI_pred\t"
            #META_SCOREPREDICTOR_BODY+=f"{polyphen2hdiv}\t{polyphen2hdiv_pred}\t{polyphen2hvar}\t{polyphen2hvar_pred}\t{primateai}\t{primateai_pred}"
            
            provean=info['PROVEAN_score'] if 'PROVEAN_score' in info else EMPTY_FIELD
            provean_pred=self.varsomeCriteria['PROVEAN_score'][4](provean)
            
            polyphen2hdiv=info['Polyphen2_HDIV_score'] if 'Polyphen2_HDIV_score' in info else EMPTY_FIELD
            polyphen2hdiv_pred=self.varsomeCriteria['Polyphen2_HDIV_score'][4](polyphen2hdiv)
            polyphen2hvar=info['Polyphen2_HVAR_score'] if 'Polyphen2_HVAR_score' in info else EMPTY_FIELD
            polyphen2hvar_pred=self.varsomeCriteria['Polyphen2_HVAR_score'][4](polyphen2hvar)
            
            primateai=info['PrimateAI_score'] if 'PrimateAI_score' in info else EMPTY_FIELD
            primateai_pred=self.varsomeCriteria['PrimateAI_score'][4](primateai)
            META_SCOREPREDICTOR_HEAD+=f"Polyphen2_HDIV\tPolyphen2_HDIV_pred\Polyphen2_HVAR\tPolyphen2_HVAR_pred\tPrimateAI\tPrimateAI_pred\t"
            META_SCOREPREDICTOR_BODY+=f"{polyphen2hdiv}\t{polyphen2hdiv_pred}\t{polyphen2hvar}\t{polyphen2hvar_pred}\t{primateai}\t{primateai_pred}"
            
            
            sift=info['SIFT_score'] if 'SIFT_score' in info else EMPTY_FIELD
            sift_pred=self.varsomeCriteria['SIFT_score'][4](sift)
            sift4g=info['SIFT4G_score'] if 'SIFT4G_score' in info else EMPTY_FIELD
            sift4g_pred=self.varsomeCriteria['SIFT4G_score'][4](sift4g)
            META_SCOREPREDICTOR_HEAD+=f"SIFT\tSIFT_pred\SIFT4G_RF\tSIFT4GF_pred\t"
            META_SCOREPREDICTOR_BODY+=f"{sift}\t{sift_pred}\t{sift4g}\t{sift4g_pred}\t"
            
            dbscsnvada=info['dbscSNV_ADA_SCORE'] if 'dbscSNV_ADA_SCORE' in info else EMPTY_FIELD
            dbscsnvada_pred=self.varsomeCriteria['dbscSNV_ADA_SCORE'][4](dbscsnvada)
            dbscsnvrf=info['dbscSNV_RF_SCORE'] if 'dbscSNV_RF_SCORE' in info else EMPTY_FIELD
            dbscsnvrf_pred=self.varsomeCriteria['dbscSNV_RF_SCORE'][4](dbscsnvrf)
            
            mcap=info['M-CAP_score'] if 'M-CAP_score' in info else EMPTY_FIELD
            mcap_pred=self.varsomeCriteria['M-CAP_score'][4](mcap)
    
            META_SCOREPREDICTOR_HEAD+=f"dbscSNV_ADA\tdbscSNV_ADA_pred\dbscSNV_RF\tdbscSNV_RF_pred\t"
            META_SCOREPREDICTOR_BODY+=f"{dbscsnvada}\t{dbscsnvada_pred}\t{dbscsnvrf}\t{dbscsnvrf_pred}\t"
            """
            
            ##############################
            ### META SCORE PREDICTOR
            ##############################
            score_pred_meta_score,score_pred_meta_path,score_pred_meta_ben,META_SCOREPREDICTOR_BODY=self.aggregateScorePredictor(info,PATHOGENIC_SCORE_LIST_METASCORE)
            """
            bayesaddaf=info['BayesDel_addAF_score'] if 'BayesDel_addAF_score' in info else EMPTY_FIELD
            bayesaddaf_pred=self.varsomeCriteria['BayesDel_addAF_score'][4](bayesaddaf)
            bayesnoaf=info['BayesDel_noAF_score'] if 'BayesDel_noAF_score' in info else EMPTY_FIELD
            bayesnoaf_pred=self.varsomeCriteria['BayesDel_noAF_score'][4](bayesnoaf)
            
            metalr=info['MetaLR_score'] if 'MetaLR_score' in info else EMPTY_FIELD
            metalr_pred=self.varsomeCriteria['MetaLR_score'][4](metalr)
            metarnn=info['MetaRNN_score'] if 'MetaRNN_score' in info else EMPTY_FIELD
            metarnn_pred=self.varsomeCriteria['MetaRNN_score'][4](metarnn)
            metasvm=info['MetaSVM_score'] if 'MetaSVM_score' in info else EMPTY_FIELD
            metasvm_pred=self.varsomeCriteria['MetaSVM_score'][4](metasvm)
            
            revel=info['REVEL_score'] if 'REVEL_score' in info else EMPTY_FIELD
            revel_pred=self.varsomeCriteria['REVEL_score'][4](revel)
            
            META_SCOREPREDICTOR_HEAD=f"BayesDel_addAF\tBayesDel_addAF_score\tBayesDel_noAF\tBayesDel_noAF_pred\tMetaLR\tMetaLR_score\tMetaRNN\tMetaRNN_score\tMetaSVM\tMetaSVM_score\tREVEL\tREVEL_score\t"
            META_SCOREPREDICTOR_BODY=f"{bayesaddaf}\t{bayesaddaf_pred}\t{bayesnoaf}\t{bayesnoaf_pred}\t{metalr}\t{metalr_pred}\t{metarnn}\t{metarnn_pred}\t{metasvm}\t{metasvm_pred}\t{revel}\t{revel_pred}\t"
            score_pred_meta=[bayesaddaf_pred,bayesnoaf_pred,metalr_pred,metarnn_pred,metasvm_pred,revel_pred]
            score_pred_meta_path=sum([1 if x>0 else 0 for x in score_pred_meta])
            score_pred_meta_ben=sum([1 if x<0 else 0 for x in score_pred_meta])
            score_pred_meta_score=sum(score_pred_meta)
            """
            ################################
            ### CONSERVATION SCORE PREDICTOR
            ################################
            score_pred_conservation_score,score_pred_conservation_path,score_pred_conservation_ben,CONSERVATION_SCOREPREDICTOR_BODY=self.aggregateScorePredictor(info, PATHOGENIC_SCORE_LIST_CONSERVATION)
            """
            phastcon=info['phastCons100way_vertebrate'] if 'phastCons100way_vertebrate' in info else EMPTY_FIELD
            phastcon_pred=self.varsomeCriteria['phastCons100way_vertebrate'][4](phastcon)
            phastCons30way_mammalian=info['phastCons30way_mammalian'] if 'phastCons30way_mammalian' in info else EMPTY_FIELD
            
            phylop=info['phyloP30way_vertebrate'] if 'phyloP30way_vertebrate' in info else EMPTY_FIELD
            phylop_pred=self.varsomeCriteria['phyloP30way_vertebrate'][4](phylop)
            phyloP30way_mammalian=info['phyloP30way_mammalian'] if 'phyloP30way_mammalian' in info else EMPTY_FIELD
            phyloP100way_vertebrate=info['phyloP100way_vertebrate'] if 'phyloP100way_vertebrate' in info else EMPTY_FIELD

            siphy=info['SiPhy_29way_logOdds'] if 'SiPhy_29way_logOdds' in info else EMPTY_FIELD
            
            GERP_NR=info['GERP++_NR'] if 'GERP++_NR' in info else EMPTY_FIELD
            GERP_RS=info['GERP++_RS'] if 'GERP++_RS' in info else EMPTY_FIELD
            
            CONS_SCOREPREDICTOR_HEAD=f"phastCons100way_vertebrate\tphastCons100way_vertebrate_pred\phyloP30way_vertebrate\phyloP30way_vertebrate_pred\tSiPhy_29way\tGERP_NR\tGERP_RS\t"
            CONS_COREPREDICTOR_BODY=f"{phastcon}\t{phastcon_pred}\t{phylop}\t{phylop_pred}\t{siphy}\t{GERP_NR}\t{GERP_RS}"
            score_pred_conservation=[phastcon_pred,phylop_pred]
            score_pred_conservation_path=sum([1 if x>0 else 0 for x in score_pred_conservation])
            score_pred_conservation_ben=sum([1 if x<0 else 0 for x in score_pred_conservation])
            score_pred_conservation_score=sum(score_pred_conservation)
            
            mvp=info['MVP_score'] if 'MVP_score' in info else EMPTY_FIELD
            mvp_pred=self.varsomeCriteria['MVP_score'][4](mvp)            
            """
            ################################
            ### OTHER SCORE PREDICTOR
            ################################
            
            #GTEx_V8_gene=info['GTEx_V8_gene'] if 'GTEx_V8_gene' in info else EMPTY_FIELD 
            GTEx_V8_tissue=info['GTEx_V8_tissue'].replace("_"," ") if 'GTEx_V8_tissue' in info else EMPTY_FIELD
            Interpro_domain=info['Interpro_domain'].replace("_"," ") if 'Interpro_domain' in info else EMPTY_FIELD
            icgc=info['ICGC_Id'] if 'ICGC_Id' in info else EMPTY_FIELD
            icgc_Occurrence=info['ICGC_Occurrence'] if 'ICGC_Occurrence' in info else EMPTY_FIELD
            
            mq=info['MQ'] if 'MQ' in info else EMPTY_FIELD
            vqslod=info['VQSLOD'] if 'VQSLOD' in info else EMPTY_FIELD
            excesshet=info['ExcessHet'] if 'ExcessHet' in info else EMPTY_FIELD
            qd=info['QD'] if 'QD' in info else EMPTY_FIELD
            
            format_field=variant.format
            DP,AD_ref,AD_alt,GT,GQ=getFORMAT_2String(format_field) if self.columnFormat_SampleNumber == 1 else getFORMAT_2String(format_field,index=0)
            if self.columnFormat_SampleNumber == 2:
                DP,AD_ref,AD_alt,GT,GQ=getFORMAT_2String(format_field, index=0)
                DP_tmp,AD_ref_tmp,AD_alt_tmp,GT_tmp,GQ_tmp=getFORMAT_2String(format_field,index=1)
                DP=(DP,DP_tmp)
                AD_ref=(AD_ref,AD_ref_tmp)
                AD_alt=(AD_alt,AD_alt_tmp)
                GT=(GT,GT_tmp)
                GQ=(GQ,GQ_tmp)
            
            #psp=list(set(self.varsomeCriteria.keys()).intersection(variant.info))
            #for psp_key in psp :
            #    value=variant.info[psp_key]
            #    variant.info["mypsp_"+psp_key]=(value,self.varsomeCriteria[psp_key][4](value))
                
            
            variant_info=variant.info['ANN']
            for trs in variant_info:
                gene=trs.gene
                transcript=trs.featureID
                isCanonical=EMPTY_FIELD if trs.isCanonical is None else trs.isCanonical 
                hgvsc=trs.HGVSc
                hgvsp=trs.HGVSp
                rank_ex=trs.rank[0] if trs.rank is not None else EMPTY_FIELD #f"{trs.rank[0]}/{trs.rank[1]}" if trs.rank is not None else EMPTY_FIELD
                rank_tt=trs.rank[1] if trs.rank is not None else EMPTY_FIELD
                annotation=trs.annotation
                type_tr=trs.transcript
                impact=trs.annotationImpact
                region=trs.region
                #HEADER_GERMLINE_STRING="gene\ttranscript\tisCanonical\tavsnp\tchr\tposition\thgvsc\thgvsp\tExone\tExT\tannotation\tannovar_FuncRefGene\tannovar_ExFunc\t"
                GERMLINE_STRING=f"{gene}\t{transcript}\t{isCanonical}\t{avsnp}\t{chromosome}\t{position}\t{hgvsc}\t{hgvsp}\t{rank_ex}\t{rank_tt}\t{annotation}\t{annovar_FuncRefGene}\t{annovar_ExFunc}\t"
                ### Genotype and QUALITY
                if self.columnFormat_SampleNumber == 1:
                    #HEADER_GERMLINE_STRING+="GT\tAD_ref\tAD_alt\tDP\tGQ\t"
                    GERMLINE_STRING+=f"{GT}\t{AD_ref}\t{AD_alt}\t{DP}\t{GQ}\t"
                elif self.columnFormat_SampleNumber == 2 :
                    #HEADER_GERMLINE_STRING+="gGT\tsGT\tgAD_ref\tsAD_ref\tgAD_alt\tsAD_alt\tgDP\tsDP\tgGQ\stGQ\t"
                    GERMLINE_STRING+=f"{GT[0]}\t{GT[1]}\t{AD_ref[0]}\t{AD_ref[1]}\t{AD_alt[0]}\t{AD_alt[1]}\t{DP[0]}\t{DP[1]}\t{GQ[0]}\t{GQ[1]}\t"
                
                #HEADER_GERMLINE_STRING+=f"gnomAD\tcadd\t{'CancerVar' if self.columnFormat_SampleNumber == 2 else 'InterVar'}\tClinvar\tClinvar_star\t"
                GERMLINE_STRING+=f"{gnomAD}\t{cadd}\t{ACMG}\t{Clinvar}\t{Clinvar_star}\t"

                ### Genotype and QUALITY
                #HEADER_GERMLINE_STRING+="omim_genecode\tomim_inheritance\tomim\tomim_pathology\tflag_lof\tlof\tflag_nmd\tnmd\ttype_tr\timpact\tregion\tannovar_NM\tref\talt\tcytoBand\t"
                GERMLINE_STRING+=f"{omim_genecode}\t{omim_inheritance}\t{omim}\t{omim_pathology}\t{flag_lof}\t{lof}\t{flag_nmd}\t{nmd}\t{type_tr}\t{impact}\t{region}\t{annovar_NM}\t{ref}\t{alt}\t{cytoBand}\t"
                
                ### Predictor: META
                #HEADER_GERMLINE_STRING+="MetaScore\tMetaScoreBen\tMetaScorePat\t"
                GERMLINE_STRING+=f"{score_pred_meta_score}\t{score_pred_meta_path}\t{score_pred_meta_ben}\t{META_SCOREPREDICTOR_BODY}"
                
                ### Predictor: INDIVIDUAL
                #HEADER_GERMLINE_STRING+="IndivScore\tIndivScoreBen\tIndivScorePat\t"
                GERMLINE_STRING+=f"{score_pred_individual_scor}\t{score_pred_individual_pat}\t{score_pred_individual_ben}\t{INDIVIDUAL_SCOREPREDICTOR_BODY}"
                
                ### Predictor: CONSERVATION
                #HEADER_GERMLINE_STRING+="ConservationScore\tConservationScoreBen\tConservationScorePat\t"
                GERMLINE_STRING+=f"{score_pred_conservation_score}\t{score_pred_conservation_path}\t{score_pred_conservation_ben}\t{CONSERVATION_SCOREPREDICTOR_BODY}"
                
                """
                HEADER_GERMLINE_STRING+="SIFT_score\tSIFT4G_score\tPolyphen2_HDIV_score\tPolyphen2_HVAR_score\tDANN\tFATHMM_score\tPROVEAN_score\tMetaSVM_score\tdbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\t"
                GERMLINE_STRING+=f"{sift}\t{sift4g}\t{polyphen2hdiv}\t{polyphen2hvar}\t{dann}\t{fathmm}\t{provean}\t{metasvm}\t{dbscsnvada}\t{dbscsnvrf}\t"
                ### Predictor: CONSERVATION
                HEADER_GERMLINE_STRING+="GERP_NR\tGERP_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphastCons100way_vertebrate\tphastCons30way_mammalian\tPrimateAI\tSiPhy\t"
                GERMLINE_STRING+=f"{GERP_NR}\t{GERP_RS}\t{phylop}\t{phyloP30way_mammalian}\t{phastcon}\t{phastCons30way_mammalian}\t{primateai}\t{siphy}\t"
                """
                #HEADER_GERMLINE_STRING+="MQ\tVSQLOD\tQD\tExcHet\tGTEx_tissue\tInterpro\ticgc\ticgc_Occurrence\tClinVar_Phenotype"
                GERMLINE_STRING+=f"{mq}\t{vqslod}\t{qd}\t{excesshet}\t{GTEx_V8_tissue}\t{Interpro_domain}\t{icgc}\t{icgc_Occurrence}\t{Clinvar_info}"
                
                self.outfile.write((",".join(self.sampleID)+"\t" if self.DUMP_SAMPLE_NAME else "")+GERMLINE_STRING+"\n")
                #self.outfile.write("{:}{:}{:}{:}{:}\t{:}\n".format(txt, newVar.replace("|", "\t") + "\t", txt_format, txt_info,flag_lof,flag_nmd))
        
        return variantsCounter

    def markSnpeffTranscript(self, variant, keepOnlyCanonical=False):
        if variant.info is not None and 'ANN' in variant.info :
            variant_info=variant.info['ANN']
            #variant_info.sort(key=lambda trs:( ({"HIGH":-1,"MODERATE":-2,"MODIFIER":-3,"LOW":-4,".":-5}[trs.annotationImpact] ), trs.cds[1] if trs.cds is not None else -1, trs.cDNA[1] if trs.cDNA is not None else -1) , reverse=True)
            filtered_list=[]
            
            if self.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT:
                for trs in variant_info :
                    if trs.gene in self._list_of_genes and trs.featureID in self._list_of_canonical_transcript :
                        trs.isCanonical=True
                        filtered_list.append(trs)
            else :
                for trs in variant_info :
                    if trs.gene in self._list_of_genes :
                        if trs.featureID in self._list_of_canonical_transcript :
                            trs.isCanonical=True
                        else :
                            trs.isCanonical=False
                        filtered_list.append(trs)
                
            variant.info['ANN']=filtered_list
        
    def keepOnlyImpactTranscripts_SnpEff(self, variant):
        variant_info=variant.info['ANN']
        filtered_list=[]
        for trs in variant_info :
            if trs.annotationImpact in self.SNPEFF_INTERESTING_IMPACT:
                filtered_list.append(trs)
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
    
    def aggregateScorePredictor(self, info, score_key_list):
        total_score=0
        score_path=0
        score_ben=0
        SCOREPREDICTOR_BODY=""
        for score_key in score_key_list :
            score_val=EMPTY_FIELD
            score_pred=EMPTY_FIELD
            if score_key in info :
                score_val=info[score_key]
                if score_key in self.varsomeCriteria :
                    score_pred=self.varsomeCriteria[score_key][4](score_val)
                    score_path+=1 if score_pred>0 else 0
                    score_ben+=1 if score_pred<0 else 0
                    total_score+=score_pred
            SCOREPREDICTOR_BODY+=f"{score_val}\t{score_pred}\t"
        return total_score,score_path,score_ben,SCOREPREDICTOR_BODY
        
    def composeTableHeader(self, line):
        ### COLLECT Information about the Sample IDs (relevant expecially for the Somatic, Trio ...)
        self.sampleID=line.strip().split('\t')[9:]
        ### ENUMERATE the SAMPLE in the OUTPUT TABLE
        self.columnFormat_SampleNumber=len(self.sampleID)
        #if self.DUMP_SAMPLE_NAME :
        #    for i in range(self.columnFormat_SampleNumber) :
        #        self.outfile.write("# Sample{:d}: {:}\n".format(i+1, self.sampleID[i]))

        """
        HEADER_GERMLINE_STRING="gene\ttranscript\tisCanonical\tavsnp\tchr\tposition\thgvsc\thgvsp\tExone\tExT\tannotation\tannovar_FuncRefGene\tannovar_ExFunc\t"
        ### Genotype and QUALITY
        HEADER_GERMLINE_STRING+="GT\tAD_ref\tAD_alt\tDP\tGQ\tgnomAD\tcadd\tACMG\tClinvar\tClinvar_star\t"
        HEADER_GERMLINE_STRING+="omim_inheritance\tomim\tomim_pathology\tflag_lof\tlof\tflag_nmd\tnmd\ttype_tr\timpact\tregion\tannovar_NM\tref\talt\tcytoBand\t"
        HEADER_GERMLINE_STRING+="SIFT_score\tSIFT4G_score\tPolyphen2_HDIV_score\tPolyphen2_HVAR_score\tDANN\tFATHMM_score\tPROVEAN_score\tMetaSVM_score\tdbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\tGERP_NR\tGERP_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphastCons100way_vertebrate\tphastCons30way_mammalian\tPrimateAI\tSiPhy\t"
        HEADER_GERMLINE_STRING+="MQ\tVSQLOD\tQD\tExcHet\tGTEx_tissue\tInterpro\ticgc\ticgc_Occurrence"
        """
        HEADER_GERMLINE_STRING="gene\ttranscript\tisCanonical\tavsnp\tchr\tposition\thgvsc\thgvsp\tExone\tExT\tannotation\tannovar_FuncRefGene\tannovar_ExFunc\t"
        ### Genotype and QUALITY
        if self.columnFormat_SampleNumber == 1:
            HEADER_GERMLINE_STRING+="GT\tAD_ref\tAD_alt\tDP\tGQ\t"
            HEADER_GERMLINE_STRING+="gnomAD\tcadd\tInterVar\tClinvar\tClinvar_star\t"
        elif self.columnFormat_SampleNumber == 2 :
            HEADER_GERMLINE_STRING+="gGT\tsGT\tgAD_ref\tsAD_ref\tgAD_alt\tsAD_alt\tgDP\tsDP\tgGQ\tsGQ\t"
            HEADER_GERMLINE_STRING+="gnomAD\tcadd\tCancerVar\tClinvar\tClinvar_star\t"

        ### Genotype and QUALITY
        HEADER_GERMLINE_STRING+="omim_genecode\tomim_inheritance\tomim\tomim_pathology\tflag_lof\tlof\tflag_nmd\tnmd\ttype_tr\timpact\tregion\tannovar_NM\tref\talt\tcytoBand\t"
        HEADER_GERMLINE_STRING+="MetaScore\tMetascorePat\tMetascoreBen\t"
        for score_key in PATHOGENIC_SCORE_LIST_METASCORE :
            HEADER_GERMLINE_STRING+=f"{score_key}\tRating\t"
        HEADER_GERMLINE_STRING+="IndivScore\tIndividualPat\tIndividualBen\t"
        for score_key in PATHOGENIC_SCORE_LIST_INDIVIDUAL :
            HEADER_GERMLINE_STRING+=f"{score_key}\tRating\t"
        HEADER_GERMLINE_STRING+="ConservationScore\tConservationPat\tConservationBen\t"
        for score_key in PATHOGENIC_SCORE_LIST_CONSERVATION :
            HEADER_GERMLINE_STRING+=f"{score_key}\tRating\t"
        #HEADER_GERMLINE_STRING+="SIFT_score\tSIFT4G_score\tPolyphen2_HDIV_score\tPolyphen2_HVAR_score\tDANN\tFATHMM_score\tPROVEAN_score\tMetaSVM_score\tdbscSNV_ADA_SCORE\tdbscSNV_RF_SCORE\t"
        #HEADER_GERMLINE_STRING+="GERP_NR\tGERP_RS\tphyloP100way_vertebrate\tphyloP30way_mammalian\tphastCons100way_vertebrate\tphastCons30way_mammalian\tPrimateAI\tSiPhy\t"
        HEADER_GERMLINE_STRING+="MQ\tVSQLOD\tQD\tExcHet\tGTEx_tissue\tInterpro\ticgc\ticgc_Occurrence\tClinVar_Phenotype"

        
        return HEADER_GERMLINE_STRING

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

        self.setAnnotation=set(self.annotationParserDictionary.keys())
        
        headertable=self.composeTableHeader(line)
        self.outfile.write(("Sample\t" if self.DUMP_SAMPLE_NAME else "")+headertable+"\n")
        
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

    def setVarsomeCriteria(self, varsomeCriteria):
        self.varsomeCriteria=varsomeCriteria
    
    def close(self):
        if not self.outfile.closed :
            self.outfile.close()

def executeTableCreation():
    parser = argparse.ArgumentParser(prog="table_new", description="Transform a VCF file annotated into a text file formatted as a table (tab separated) extracting a subset of fields. ")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="Input annotated VCF data file (produced by SPNEff or VEP).", required=True, default=None)
    parser.add_argument('-s','--dump_sample',action='store_true', help="Dump the name of the sample [DEFAULT: False].", default=False)
    parser.add_argument('-only','--only_transcript_with_impact',action='store_true', help="Flag excluding all transcript that do no have an impact HIGH or MODERATE [DEFAULT: False].", default=False)
    parser.add_argument('-canonical','--keep_only_canonical_transcript',action='store_true', help="Flag excluding all non canonical transcript (is canonical the transcript with the larger CDS or for non coding the larger cDNA) [DEFAULT: False].", default=False)
    
    ### Read CSV - Interesting Genes Names
    parser.add_argument('-gl','--gene_list',action='store',type=str, help="Path of the file containing the genes od interest [DEFAULT: None].", default=None)
    parser.add_argument('-gls','--gene_list_separator',action='store',type=str,help="Column separator in the gene_list file [DEFAULT: '\t'].", required=False, default="\t")
    parser.add_argument('-glc','--gene_list_column',action='store',type=int, help=f"Index of the column containing the GENE names convention is to start column count from zero [DEFAULT: 0].", default=0)
    parser.add_argument('-tlc','--gene_list_transcript_column',action='store',type=int, help=f"Index of the column containing the TRANSCRIPT names convention is to start column count from zero [DEFAULT: 1].", default=1)
    parser.add_argument('-gln','--gene_list_noheader',action='store_true', help="GeneList file contains header [DEFAULT: False].", default=False)

    parser.add_argument('-vs','--varsome',action='store',type=str, help=f"Table of criteria for the pathogenic score predictor [DEFAULT: {DEFAULT_VARSOME_CRITERIA_TABLE}].", default=DEFAULT_VARSOME_CRITERIA_TABLE)
    
    ### PARAMETERs
    parser.add_argument('-o','--output',action='store',type=str,help="Output file path used to store the tabular results .", required=False, default="Table_output.tsv")
    arg=parser.parse_args()

    INPUT=arg.input
    ### Load the field for INFO
    output=arg.output
    
    varsomeCriteria=loadVarsomeCriteriaTable(arg.varsome)
    print("INPUT FILE: {:}".format(INPUT))
    reader=readerVCF_AnnotationFilter(INPUT, output=output)
    reader.DUMP_SAMPLE_NAME=arg.dump_sample
    reader.KEEP_ONLY_TRANSCRIPT_WITH_IMPACT=arg.only_transcript_with_impact
    reader.KEEP_ONLY_SNPEFF_CANONICAL_TRANSCRIPT=arg.keep_only_canonical_transcript
    
    ### Load the GENE NAMES to be KEPT in the TABLE
    geneNameList=None
    if arg.gene_list is not None and arg.gene_list != "" :
        geneNameList=loadCSV_Column(arg.gene_list, arg.gene_list_separator, column_key=arg.gene_list_column, column_value=arg.gene_list_transcript_column, header=(not arg.gene_list_noheader))
    else :
        print("No specific gene list has been loaded all canonical genes will be used.")
    reader.setGeneNames(set(geneNameList.keys()))
    reader.setCanonicalTranscript(set(geneNameList.values()))
    reader.setVarsomeCriteria(varsomeCriteria)
    
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


################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    executeTableCreation()