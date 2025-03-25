import csv, re,sys

#0. read the user specified evidence file
def loadUserEvidenceDB(db_path):
    user_evidence_dict={}
    try:
        fh=open(db_path, "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            print("EvidenceFile:\t"+line2)
            cls2=line2.split('\t')
            if len(cls2)>1:
                keys=cls2[0]+"_"+cls2[1]+"_"+cls2[2]+"_"+cls2[3]+"_"+cls2[4]
                keys=re.sub("[Cc][Hh][Rr]","",keys)
                user_evidence_dict[keys]=cls2[5].upper()
                print("EvidenceDict:\t"+user_evidence_dict)
                print("%s %s\n" %(keys,user_evidence_dict[keys]))
    except IOError:
        print("Error: can\'t read the user specified evidence file %s" % db_path)
    else:
        fh.close()
    return user_evidence_dict


#1.LOF gene list
def loadLofDB(db_path):
    lof_genes_dict={}
    try:
        fh = open(db_path, "r")
        strs = fh.read()
        count=0
        for line2 in strs.split('\n'):
            if len(line2)>0 and line2[0]!='#' :
                count+=1
                cls2=line2.split('\t')
                if len(cls2[0])>1:
                    lof_genes_dict[cls2[0]]='1'
    except IOError as err:
        print(f"ERROR!!!!! {err}")
        print("Error: can\'t read the LOF genes file %s" % db_path)
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()
    return lof_genes_dict

#1. cancervar_markers
def loadCancervarDB(db_path):
    cancervar_d=[]
    cancervar_markers_dict={}
    try:
        with open(db_path) as fh:
            reader = csv.reader(filter(lambda row: row[0]!='#', fh), delimiter="\t", )
            cancervar_d = list(reader)
    
    except IOError:
        print("Error: can\'t read the cancervar_markers file %s" % db_path)
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()
        #print len(cancervar_d)
        for ii,row in enumerate(cancervar_d):
            gene=row[0]
            mut_type=row[1]
            mut=row[2]
            mut_types=mut_type.split(';');  # BIA CNA EXPR FUS MUT
            mut_alts=mut.split(';');
            for jj,mut_type_el in enumerate(mut_types):
                mut_list=mut_alts[jj].split(':')
                if(len(mut_list)>1):
                    gene1=mut_list[0]
                    mutb=mut_list[1]
                else:
                    gene1=gene
                    mutb=mut_alts[jj]
                if(mut_types[jj]=="FUS"):                    
                    mutt="fus_"+re.sub('__','-',mut_alts[jj])
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="EXP" or mut_types[jj]=="EXPR"):                    
                    mutt=gene1+"_"+"expr_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="CNV" or mut_types[jj]=="CNA"):                    
                    mutt=gene1+"_"+"cna_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="BIA"):                    
                    mutt=gene1+"_"+"bia_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                elif(mut_types[jj]=="OTH"):                    
                    mutt=gene1+"_"+"oth_"+mutb
                    key=mutt
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                #start from mutb
                elif(mut_types[jj]=="MUT"):                    
                    if(mutb=="any mutation"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any insertion"): 
                        mutt="insertion_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any deletion"): 
                        mutt="deletion_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any frameshift" or mutb=="frameshift"): 
                        mutt="frameshift_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any indel"): 
                        mutt="indel_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any missense"): 
                        mutt="missense_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="mutant"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="mutation" or mutb=="ALTERATION"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="MUTATION"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="oncogenic mutation"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="Oncogenic Mutations"): 
                        mutt="mut_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    elif(mutb=="any nonsense"): 
                        mutt="nonsense_any"
                        key=gene+'_'+mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
    
                    # format as " exon(s) 10, 20, 21 any"
                    elif(  re.findall('exon\(s\) ', mutb, flags=re.IGNORECASE)):
                        str2=mutb.split('exon(s) ',1)[1];
                        if( not  re.findall(',', mutb, flags=re.IGNORECASE)): # for single position
                            poss,mutc=str2.split(' ',1);
                            #mut0=re.sub('"','',mutb)
                            poss_t=int(poss)
                            muts="exon_"+str(poss_t)+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                        else: # for multiple positions
                            list_pos=str2.split(',');
                            list_pos_size=len(list_pos)
    
                            list_tt=list_pos[list_pos_size-1].split(' ')
                            list_tt_size=len(list_tt)
    
                            mutc=list_tt[list_tt_size-1]
                            pos_end=str(int(list_tt[list_tt_size-2]))
                            if(mutc=="mutation"): mutc="any"
                            muts="exon_"+pos_end+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            for jj in range(0,list_pos_size-1):
                                pos_be=str(int(list_pos[jj]))
                                if(mutc=="mutation"): mutc="any"
                                muts="exon_"+pos_be+"_"+mutc
                                key=gene+'_'+muts
                                default_s=''
                                cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                                #print key,ii,cancervar_markers_dict[key]
                    # " codon(s) 289, 596, 598 any"   codon(s) 132 any
                    elif(  re.findall('codon\(s\) ', mutb, flags=re.IGNORECASE)):
                        str2=mutb.split('codon(s) ',1)[1];
                        if( not  re.findall(',', mutb, flags=re.IGNORECASE)): # for single position
                            poss,mutc=str2.split(' ',1);
                            #mut0=re.sub('"','',mutb)
                            poss_t=int(poss)
                            muts="codon_"+str(poss_t)+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                        else: # for multiple positions
                            list_pos=str2.split(',');
                            list_pos_size=len(list_pos)
    
                            list_tt=list_pos[list_pos_size-1].split(' ')
                            list_tt_size=len(list_tt)
    
                            mutc=list_tt[list_tt_size-1]
                            pos_end=str(int(list_tt[list_tt_size-2]))
                            muts="codon_"+pos_end+"_"+mutc
                            key=gene+'_'+muts
                            default_s=''
                            cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                            for jj in range(0,list_pos_size-1):
                                pos_be=str(int(list_pos[jj]))
                                muts="codon_"+pos_be+"_"+mutc
                                key=gene+'_'+muts
                                default_s=''
                                cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    else:
                        mutt=gene1+"_"+mutb
                        key=mutt
                        default_s=''
                        cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                else:
                    key=gene1+"_"+mutb
                    default_s=''
                    cancervar_markers_dict[key]=str(ii)+","+cancervar_markers_dict.get(key,default_s)
                    
    return cancervar_d,cancervar_markers_dict

#10 cancer_pathway=%(database_cancervar)s/cancer_pathway.list
def loadCancerPathways(db_path):
    cancer_pathway_dict={}
    try:
        fh = open(db_path, "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            if len(line2)>0 and line2[0]!='#' :
                cls2=line2.split('\t')
                if len(cls2[0])>1:
                    #cancer_pathway_dict[cls2[0]]='1'
                    cancer_pathway_dict[cls2[0]]=cls2[1]
    except IOError:
        print("Error: can\'t read the cancer_pathway genes file %s" % db_path)
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()
    return cancer_pathway_dict

#11 cancers_genes=%(database_cancervar)s/cancers_genes.list
def loadCancerGenes(db_path):
    cancers_gene_dict={}
    try:
        fh = open(db_path, "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            if len(line2)>0 and line2[0]!='#' :
                cls2=line2.split('\t')
                if len(cls2[0])>1:
                    cancers_gene_dict[cls2[0]]='1'
    except IOError:
        print("Error: can\'t read the cancers diseases genes file %s" % db_path)
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()
    return cancers_gene_dict

#12 cancers_types=%(database_cancervar)s/cancervar.cancer.types
def loadCancerTypes(db_path):
    cancers_types_dict={}
    try:
        fh = open(db_path, "r")
        strs = fh.read()
        for line2 in strs.split('\n'):
            if len(line2)>0 and line2[0]!='#' :
                cls2=line2.split('/')
                if len(cls2)>1:
                    cancers_types_dict[cls2[1]]=cls2[0];
    except IOError:
        print("Error: can\'t read the cancers types file %s" % db_path)
        print("Error: Please download it from the source website")
        sys.exit()
        return
    else:
        fh.close()
    return cancers_types_dict

