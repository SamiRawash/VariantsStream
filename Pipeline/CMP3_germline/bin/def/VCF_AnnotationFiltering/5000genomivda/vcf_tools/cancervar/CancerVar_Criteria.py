import re
from vcf_tools.cancervar.CancerVar_LoadDB import loadCancerGenes,loadCancerPathways,loadCancerTypes,loadCancervarDB,loadLofDB,loadUserEvidenceDB

GNOMAD_FIELDS=set(["AF","AF_raw","AF_male","AF_female","AF_afr","AF_ami","AF_amr","AF_asj","AF_eas","AF_fin","AF_nfe","AF_oth","AF_sas"])
GNOMAD_FIELDS_RARE_TRASHOLD=0.01
GNOMAD_FIELDS_MAF_CUTOFF=0.0001

SCORE_PREDICTOR_dbscSNV_CUTOFF=0.6
SCORE_PREDICTOR_SIFT_CUTOFF=0.05
SCORE_PREDICTOR_METASVM_CUTOFF=0
SCORE_PREDICTOR_GERP_CUTOFF=2

### TIER - CBP TRESHOLD
##CANCERVAR_BPS=["Pathogenic","Likely_pathogenic","Benign/Likely_benign","Uncertain_significance"]
CANCERVAR_BPS=["Tier_I_strong","Tier_II_potential","Tier_IV_benign","Tier_III_Uncertain"]
CBP_TRESHOLD_BENIGN_LIKELYBENIGN_MAX=2
CBP_TRESHOLD_LIKELY_PATHOGENIC_MIN=8
CBP_TRESHOLD_PATHOGENIC_MIN=11

cancervar_d=None
cancervar_markers_dict=None
lof_genes_dict=None
cancer_pathway_dict=None
cancers_gene_dict=None
cancer_types_dict=None

user_evidence_dict={}

def initialize_all_database(marker_db, lof_db, pathways_db, cancer_gene_db, cancer_type_db, user_evidence_db=None):
    global cancervar_d
    global cancervar_markers_dict
    cancervar_d,cancervar_markers_dict=loadCancervarDB(marker_db) #if os.path.isfile(marker_db) else {}
    
    global lof_genes_dict
    lof_genes_dict=loadLofDB(lof_db)
    
    global cancer_pathway_dict
    cancer_pathway_dict=loadCancerPathways(pathways_db)
    
    global cancers_gene_dict
    cancers_gene_dict=loadCancerGenes(cancer_gene_db)
    
    global cancer_types_dict
    cancer_types_dict=loadCancerTypes(cancer_type_db)
    
    global user_evidence_dict
    user_evidence_dict={} if user_evidence_db is None else loadUserEvidenceDB(user_evidence_db)


def searchReferenceInCancervarDB(gene_tr, func, exonfunc, AAChange, cancer_type):
    #line_tmp=f"{gene_tr} {func} {exonfunc} {cancer_type}"
    #line_tmp2=cls[Funcanno_flgs['AAChange.refGene']]
    line_tmp2_sp=re.split(';|,',AAChange)
    exon_func=exonfunc
    marker_key0=gene_tr+"_"+"mut_any"
    out_list=[]
    for cls0 in line_tmp2_sp:
        cls0_1=cls0.split(':')
        if(len(cls0)>1 and len(line_tmp2_sp)>0 and len(cls0_1)>=4 ):
            cls0_1=cls0.split(':')
            gene=cls0_1[0]
            #transc=cls0_1[1]
            exon=re.sub('exon','',cls0_1[2])
            #cdna_change=re.sub('c.','',cls0_1[3])
            amino_change="NA";
            if(len(cls0_1)==5):amino_change=re.sub('p.','',cls0_1[4])
            ltt=len(amino_change)
            codon_num="NA";
            if(len(cls0_1)==5):codon_num=amino_change[1:ltt-1]
            #print gene, transc,exon,cdna_change,amino_change,cancer_type 
            marker_key=gene+"_"+amino_change
            marker_key0=gene+"_"+"mut_any"
            marker_key1=gene+"_"+"exon_"+exon+"_any"
            marker_key2=gene+"_"+"codon_"+codon_num+"_any"
            marker_key00=""       
            marker_key11=""       
            marker_key22=""       
            if(AAChange=="nonsynonymous SNV"):
                marker_key00=gene+"_"+"missense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_missense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_missense"
                
            if exon_func.find("frameshift")>=0 and exon_func.find("nonframe")<0 : 
                marker_key00=gene+"_"+"frameshift_any"
                marker_key11=gene+"_"+"exon_"+exon+"_frameshift"
                marker_key22=gene+"_"+"codon_"+codon_num+"_frameshift"

            if exon_func.find("stopgain")>=0 or exon_func.find("stoploss")>=0:  
                marker_key00=gene+"_"+"nonsense_any"
                marker_key11=gene+"_"+"exon_"+exon+"_nonsense"
                marker_key22=gene+"_"+"codon_"+codon_num+"_nonsense"

            if exon_func.find("deletion")>=0 :  
                marker_key00=gene+"_"+"deletion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_deletion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_deletion"
            if exon_func.find("insertion")>=0 :  
                marker_key00=gene+"_"+"insertion_any"
                marker_key11=gene+"_"+"exon_"+exon+"_insertion"
                marker_key22=gene+"_"+"codon_"+codon_num+"_insertion"
            # caution: frameshift insertion or deletion 
            #deletion frameshift indel missense nonsense_any
            # nonsynonymous SNV;stopgain/stoploss(nonsense);frameshift insertion/deletion/substitution; nonframeshift insertion/deletion/substitution;
            default_s=""
            add_list=cancervar_markers_dict.get(marker_key,default_s)+cancervar_markers_dict.get(marker_key0,default_s)+cancervar_markers_dict.get(marker_key00,default_s)+cancervar_markers_dict.get(marker_key1,default_s)+cancervar_markers_dict.get(marker_key11,default_s)+cancervar_markers_dict.get(marker_key2,default_s)+cancervar_markers_dict.get(marker_key22,default_s)
            out_list.append(add_list)
    return out_list
    
    
def check_Thera(cancer_type,related_literature):
    '''Therapeutic: 
    2 FDA approved or investigational with strong evidence
    1 FDA approved for different tumor type; investigational therapies with some evidence
    0 Cancer genes: none
    0 None  
    # the function of check_Thera  check_Diagno check_Progno are similar, should combine as one function,
      but in future for  different purpose and criteria, leave three functions.
    '''
    level=0 # ABCD
    out_list=""
    #searchReferenceInCancervarDB
    for cgi_list in related_literature :
        level=0 # ABCD
        #out_list=""
        for i in cgi_list.split(","):
            #print i
            if(len(i)>0):
                pos=int(i)
                if(cancervar_d[pos][9]=="Therapeutic"):
                    out_list=str(pos)+","+out_list
                    t_level=cancervar_d[pos][10]
                    #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                    if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                        if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE) ):
                            level=1;
                        elif level<2 : level=2 ;
                    elif level < 1 and ( t_level=='C' or t_level=='D' ) : level=1;
                    #print cancer_type,cancervar_d[pos][8]
                #print gene,level,"Therapeutic"
    Thera=level 
    if(out_list==""):  out_list="."

    out_list_t=list(set(out_list.split(","))) 
    str_out_t = ",".join(out_list_t)
    out_list=str_out_t;
    return(Thera,out_list)
    #print line_tmp

def check_Diagno(cancer_type,related_literature):
    '''Diagnostic: 
    In PG or reported evidence with consensus
    not in PG but with convincing published data
    Cancer genes: none
    None
    ''' 
    level=0 # ABCD
    out_list=""
    for cgi_list in related_literature :
        level=0 # ABCD
        #out_list=""
        for i in cgi_list.split(","):
            #print i
            if(len(i)>0):
                pos=int(i)
                if(cancervar_d[pos][9]=="Diagnostic"):
                    t_level=cancervar_d[pos][10]
                    out_list=str(pos)+","+out_list
                    #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                    if(t_level=='A'  and level < 2 ): level=2;
                    if(t_level=='B'  and level < 2 ): level=2;
                    if(t_level=='C'  and level < 1 ): level=1;
                    if(t_level=='D'  and level < 1 ): level=1;
                    if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                        if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                            level=1;
        
        #print gene,level
    Diagno=level; 
    if(out_list==""):  out_list="."

    out_list_t=list(set(out_list.split(",")))
    str_out_t = ",".join(out_list_t)
    out_list=str_out_t;

    return(Diagno,out_list)

def check_Progno(cancer_type,related_literature):
    '''Prognostic: 
    In PG or reported evidence with consensus
    not in PG but with convincing published data
    Cancer genes: none
    None
    ''' 
    level=0
    out_list=""
    for cgi_list in related_literature :
        level=0 # ABCD
        #out_list=""
        for i in cgi_list.split(","):
            #print i
            if(len(i)>0):
                pos=int(i)
                if(cancervar_d[pos][9]=="Prognostic"):
                    t_level=cancervar_d[pos][10]
                    out_list=str(pos)+","+out_list
                    #print gene,pos,t_level,cancervar_d[pos][10],cancervar_d[pos][9]
                    if(t_level=='A'  and level < 2 ): level=2;
                    if(t_level=='B'  and level < 2 ): level=2;
                    if(t_level=='C'  and level < 1 ): level=1;
                    if(t_level=='D'  and level < 1 ): level=1;
                    if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                        if(t_level=='A' or t_level=='B'):  # did not find the specific type, decrease level
                            level=1;
            #print gene,level

    Progno=level
    if(out_list==""):  out_list="."
    
    out_list_t=list(set(out_list.split(",")))
    str_out_t = ",".join(out_list_t)
    out_list=str_out_t;

    return(Progno,out_list)


def check_Mut(gene, line_tmp, dbscSNV_RF_SCORE, dbscSNV_ADA_SCORE):
    '''Mutation type:
    1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
    1 Activating, LOF (missense, nonsense, indel, splicing), CNVs, fusions
    0 Functionally unknown; mostly missense, in-frame indels; less commonly,other types
    0 Functionally benign or unknown; mostly missense; less commonly, other types
    '''

    #cls=line.split('\t')
    funcs_tmp=["nonsynonymous","missense","nonsense","frameshift","splic","stopgain","stoplost","CNV","fusion"]
    funcs_tmp2="nonframe"
    Mut=0 # 0 for Tire3/4; 1 for Tire 1/2
    VS_t1=0
    VS_t2=0
    VS_t3=0
    #dbscSNV_cutoff=0.6    #either score(ada and rf) >0.6 as splicealtering
    # Funcanno_flgs={'Func.refGene':0,'ExonicFunc.refGene':0
    for fc in funcs_tmp:
        if line_tmp.find(fc)>=0 and line_tmp.find(funcs_tmp2)<0 :
            VS_t1=1
            break
    try:
        if lof_genes_dict[ gene ] == '1' :
            VS_t2=1
    except KeyError:
        VS_t2=0
    else:
        pass
    # begin check the site in  affect the splicing
    try:
        if float(dbscSNV_RF_SCORE)>SCORE_PREDICTOR_dbscSNV_CUTOFF or float(dbscSNV_ADA_SCORE)>SCORE_PREDICTOR_dbscSNV_CUTOFF:
            VS_t3=1
    except ValueError:
        pass
    else:
        pass
    if VS_t1 !=0 and VS_t2 != 0 :
        Mut=1
    if VS_t3 !=0 and VS_t2 != 0:
        Mut=1
    #print "mut=",Mut
    return (Mut)

def check_VF(variant):
    '''Variant frequencies
    1 Mostly mosaic
    1 Mostly mosaic
    0 Mosaic or nonmosaic
    0 Mostly nonmosaic (VAF, approximately 50% or 100%)
    '''
    VF=0;
    """
    try:
        if 'AD' in variant.format :
            ad=variant.format['AD']
            if isinstance(ad, tuple):
                grm_ad_ref=
                grm_ad_alt=
                som_ad_ref=
                som_ad_alt=
            else :
                ad_ref=
                ad_alt=
                
    except Exception :
        pass
    """
    return(VF)

def check_PotG(variant):
    '''Potential germline
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF approximately 50% or 100%)
    Mostly nonmosaic (VAF, approximately 50% or 100%)
    '''
    PotG=0
    """
    try:
        if 'AD' in variant.format :
            ad=variant.format['AD']
            if isinstance(ad, tuple):
                grm_ad_ref=
                grm_ad_alt=
                som_ad_ref=
                som_ad_alt=
            else :
                ad_ref=
                ad_alt=
                
    except Exception :
        pass
    """
    return (PotG)

def check_PopD(dict_freq):
    '''Population database: ESP, dbSNP, 1000Genome, ExAC, gnomad
    1  Absent or extremely low MAF
    1  Absent or extremely low MAF
    1  Absent or extremely low MAF
    0  MAF>0.01% in the general population; or high MAF in some ethnic populations
    '''
    PopD=0;
    #Freqs_3pops={'1000g2015aug_all':0,'esp6500siv2_all':0,'ExAC_ALL':0,'AF_nfe':0}
    Freqs_3pops={'AF':0,'AF_nfe':0}
    
    tt=1;
    for val in dict_freq.values():
        if(val!='.'): # absent in all  controls
            tt=tt*0;
    if tt==1:
        PopD=1

    for key in Freqs_3pops.keys():
        try:
            if (key in dict_freq and dict_freq[key]!='.') :
                if float(dict_freq[key])>GNOMAD_FIELDS_RARE_TRASHOLD : PopD=0 #  MAF>1%
                elif float(dict_freq[key])<GNOMAD_FIELDS_MAF_CUTOFF : PopD=1  #  extremely low MAF
        except ValueError:
            pass
        else:
            pass

    return (PopD)

def check_GermD(clinvar):
    '''Germline database: HGMD, ClinVar
    2 May or may not be present
    2 May or may not be present
    0 Absent or downgraded from pathogenic to VUS
    1 Absent or present but downgraded to benign/likely benign
    '''
    GermD=0
    if clinvar is None or clinvar == "." : GermD=0
    elif clinvar.find("enign")<0 and clinvar.find("athogenic")>=0:
        GermD=1
    elif clinvar.find("ikely benign")>=0 or clinvar.find("enign")>=0:
        GermD=-1
    elif clinvar.find("enign")>=0 and clinvar.find("athogenic")>=0:
        GermD=0
    if clinvar.find("ncertain significance")>=0 :
        GermD=0
    #print "GermD=",GermD,cls[Funcanno_flgs['CLNSIG']]
    return(GermD)

def check_SomD(cosmic,icgc):
    '''Somatic database: COSMIC, My Cancer Genome, TCGA
    2 Most likely present
    1 Likely present
    0 Absent or present without association to specific tumors (potential germline VUS); present but in very few cases
    0 Absent or present without association to specific tumors (potential rare germline polymorphism)
    ''' # cosmic91    ID=COSM12560;OCCURENCE=60(haematopoietic_and_lymphoid_tissue)
    #ICGC_Id ICGC_Occurrence MU31370893  COCA-CN|1|187|0.00535,PRAD-CA|1|124|0.00806,SKCA-BR|1|66|0.01515,MELA-AU|2|183|0.01093
    SomD=0;
    if cosmic!="." and icgc!=".":
        SomD=2 
    elif cosmic!="." or icgc!=".":
        SomD=1 
    #elif cosmic=="." and icgc==".":
    #    SomD=0
    return(SomD)
          
def check_PreP(sift,polyphen2,fathmm,mtass,gerp):
    '''Predictive software: SIFT, PolyPhen2, MutTaster, CADD, MetaSVM  GERP++, 
    2 Mostly damaging; information to be used for reference only >6
    1 Mostly damaging; information to be used for reference only >3
    0 Variable; information to be used for reference only
    -1 Mostly benign; information to be used for reference only
    '''
    # MetaSVM SIFT Polyphen2_HDIV MetaLR FATHMM  MutationAssessor
    # remove MetaSVM and MetaLR
    dam=0;
    var=0;
    ben=0;
    PreP=0;

    try:
        if float(sift) >= SCORE_PREDICTOR_SIFT_CUTOFF:
            ben=ben+1
        else:
            dam=dam+1
    except ValueError:  # the sift absent means many:  synonymous indel  stop, but synonymous also is no impact
        var=var+1
    else:
        pass


    if polyphen2 == "P" or polyphen2 == "D":
        dam=dam+1
    elif polyphen2 == "B" :
        ben=ben+1
    elif polyphen2 == "." :
        var=var+1

    if fathmm == "D":
        dam=dam+1
    elif fathmm == "T" :
        ben=ben+1
    elif fathmm == "." :
        var=var+1

    if mtass == "H" or mtass == "M":
        dam=dam+1
    elif mtass == "L" or mtass == "N":
        ben=ben+1
    elif mtass == "." :
        var=var+1
    
    if gerp == ".": 
        var=var+1
    else:
        if float(gerp)>= SCORE_PREDICTOR_GERP_CUTOFF:
            dam=dam+1
        else:
            ben=ben+1

    if dam==ben or var>2: PreP=0;
    elif dam >4: PreP=2;
    elif dam >2: PreP=1;
    elif ben >2: PreP=-1;

    return(PreP)

def check_Path(gene):
    '''Pathway involvement
    2 Disease-associated pathways
    1 Involve disease-associated pathways or pathogenic pathways
    0 May or may not involve disease-associated pathways
    0 May or may not involve disease-associated pathways
    '''
    Path=0;
    
    if cancer_pathway_dict is None :
        return Path

    try:
        if cancer_pathway_dict[gene] != '' :
            Path=1
    except KeyError:
        pass
    else:
        pass
    try:
        if cancers_gene_dict[gene] == '1' :
            Path=1
    except KeyError:
        pass
    else:
        pass
    #print "Path=",Path
    return(Path)

def check_Pubs(cancer_type,related_literature):
    '''Publications: functional study, population study, other
    Therapeutic/Diagnostic/Prognostic: reported evidence with consensus
    Therapeutic: evidence of using FDA-approved therapies for different tumor types; phase 2 or 3 clinical trials for investigational therapies; Diagnostic/Prognostic: multiple lines of reported evidence without consensus
    None or no convincing evidence to determine clinical/biological significance
    Reported evidence supportive of benign/likely benign; or none
    '''
    level=0
    for cgi_list in related_literature :
        #print gene,add_list
        level=0 # ABCD
        for i in cgi_list.split(","):
            #print i
            if(len(i)>0):
                pos=int(i)
                if(cancervar_d[pos][7]!=";" and cancervar_d[pos][7]!=""):
                    level=2
                    if( not  re.findall(cancer_type, cancervar_d[pos][8], flags=re.IGNORECASE)   ):
                        level=1;  # did not find the specific type, decrease level
            #print gene,level

    Pubs=level
    return(Pubs)

def classfy(CBP,key):
    # CBP[0]:Therapeutic(2100) CBP[1]:Diagno  CBP[2]:Progno CBP[3]:Mutation(1100) CBP[4]:Variant_freq CBP[5]:Potential_germ
    # CBP[6]: Populatio(1100)  CBP[7]:Germline dat(2201) CBP[8]:Somatic dat(2100) CBP[9]:Predict_dama(2201) 
    # CBP[10]:  Path(2100)  CBP[11] : Pubs

    #begin process the user's flexible grade  to get the final interpretation
    if user_evidence_dict is not None and len(user_evidence_dict)<1:
        try:
            evds=user_evidence_dict[key] #PS1=1;PM1=1;BA1=1;PVS1 PP BS BP
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and  re.findall('grade', evd_t[0], flags=re.IGNORECASE) ):
                    #10  104353782   G   A   PVS1=1;PP1=1;PM3=1;grade_PP1=2;
                    if int(evd_t[1])<=3:
                        if(evd_t[0].find('CBP')!=-1):
                            t=evd_t[0].find('CBP');
                            tt=evd_t[0];
                            tt3=int(tt[t+3:t+4])
                            if(t<len(evd_t[0])-2 and tt3<=12 ): CBP[tt3-1]=int(evd_t[1])
        except KeyError:
            pass
        else:
            pass
    
    BPS_out=3 # BPS=[3]:Uncertain significance
    CBP_sum=sum(CBP)
    if(CBP_sum<=CBP_TRESHOLD_BENIGN_LIKELYBENIGN_MAX):  BPS_out=2 # benign
    elif(CBP_sum>=CBP_TRESHOLD_LIKELY_PATHOGENIC_MIN and CBP_sum<CBP_TRESHOLD_PATHOGENIC_MIN):  BPS_out=1 # potential path
    elif(CBP_sum>=CBP_TRESHOLD_PATHOGENIC_MIN):  BPS_out=0 # strong path

    return CBP_sum,CANCERVAR_BPS[BPS_out]

def assign(variant, cancer_type, evidence_db=None):
    CBP=[0,0,0,0,0,0,0,0,0,0,0,0]
    #chrom=chromosomeCode2sting(variant.chromosome) #re.sub("[Cc][Hh][Rr]","",variant.chromosome)    
    position_start=variant.position
    position_end=variant.position
    gene=variant.info['Gene.refGene_latest'] if 'Gene.refGene_latest' in variant.info else None
    func=variant.info['Func.refGene_latest'] if 'Func.refGene_latest' in variant.info else None
    exonfunc=variant.info['ExonicFunc.refGene_latest'] if 'ExonicFunc.refGene_latest' in variant.info else None
    AAChange=variant.info['AAChange.refGene_latest'] if 'AAChange.refGene_latest' in variant.info else None
    line_tmp=f"{func} {exonfunc}"
    
    dbscSNV_RF_SCORE=variant.info['dbscSNV_RF_SCORE'] if 'dbscSNV_RF_SCORE' in variant.info else -1
    dbscSNV_ADA_SCORE=variant.info['dbscSNV_ADA_SCORE'] if 'dbscSNV_ADA_SCORE' in variant.info else -1
    dict_freq={key: variant.info[key] for key in GNOMAD_FIELDS.intersection(variant.info)}
    clinvar=variant.info['CLNSIG'] if 'CLNSIG' in variant.info else None
    cosmic=variant.info['cosmic91'] if 'cosmic91' in variant.info else None
    icgc=variant.info['ICGC_Id'] if 'ICGC_Id' in variant.info else None
    
    sift=variant.info['SIFT_score']  if 'SIFT_score' in variant.info else None
    polyphen2=variant.info['Polyphen2_HDIV_pred']  if 'Polyphen2_HDIV_pred' in variant.info else None
    fathmm=variant.info['FATHMM_pred']  if 'FATHMM_pred' in variant.info else None
    mtass=variant.info['MutationAssessor_pred']  if 'MutationAssessor_pred' in variant.info else None
    gerp=variant.info['GERP++_RS']  if 'GERP++_RS' in variant.info else None
    
    related_literature=searchReferenceInCancervarDB(gene,func,exonfunc,AAChange, cancer_type)
    Therapeutic=check_Thera(cancer_type,related_literature) #(gene,func,exonfunc,AAChange, cancer_type)
    CBP[0],Therap_list=Therapeutic

    Diagnosis=check_Diagno(cancer_type,related_literature) #(gene,func,exonfunc,AAChange, cancer_type)
    CBP[1],Diag_list=Diagnosis

    Prognosis=check_Progno(cancer_type,related_literature) #(gene,func,exonfunc,AAChange, cancer_type)
    CBP[2],Prog_list=Prognosis
    
    ### Mutation Type
    CBP[3]=check_Mut(gene, line_tmp, dbscSNV_RF_SCORE, dbscSNV_ADA_SCORE)
    ### Variant_freq
    CBP[4]=check_VF(variant)
    ### Potential_germ
    CBP[5]=check_PotG(variant)
    ### Population_data
    CBP[6]=check_PopD(dict_freq)
    ### Germline_data
    CBP[7]=check_GermD(clinvar)
    ### Somatic_data
    CBP[8]=check_SomD(cosmic,icgc)
    ### Predict_pathoge
    CBP[9]=check_PreP(sift,polyphen2,fathmm,mtass,gerp)
    ### Pathway_invol
    CBP[10]=check_Path(gene)
    ### Publications
    CBP[11]=check_Pubs(cancer_type,related_literature) #(gene,func,exonfunc,AAChange, cancer_type)

    ### BEGING process the USER's EVIDENCE file
    key=f"{variant.chromosome}_{position_start}_{position_end}_{variant.reference}_{variant.alteration}"
    #cls[Allels_flgs['Chr']]+"_"+cls[Allels_flgs['Start']]+"_"+cls[Allels_flgs['End']]+"_"+cls[Allels_flgs['Ref']]+"_"+cls[Allels_flgs['Alt']]
    if evidence_db is not None:
        try:
            evds=user_evidence_dict[key] #CBP1=1;CBP2=1;Cancer=LUAD;
            for evd in evds.split(';'):
                evd_t=evd.split('=')
                if(len(evd_t)>1 and (not re.findall('grade', evd_t[0], flags=re.IGNORECASE)) ):
                    if int(evd_t[1])<=1:
                        #print ("%s %s %s " %(key,evd_t[1],evd_t[0]))
                        if(evd_t[0].find('CBP')!=-1):
                            t=evd_t[0].find('CBP');
                            tt=evd_t[0];
                            tt3=int(tt[t+3:t+4])
                            if(t<len(evd_t[0])-2 and tt3<=12 ): CBP[tt3-1]=int(evd_t[1])
        except KeyError:
            pass
        else:
            pass
    ### END process the USER's EVIDENCE file

    score,tier=classfy(CBP,key)

    #BP=BP_out
    return(tier,score,CBP,Therap_list,Diag_list,Prog_list)
    #return(tier,score,CBP)
