from abc import ABC, abstractmethod

from vcf_tools.TypeAnnotation import typeAnnotation, typeAnnotation_LOF  # , typeAnnotation_GnomAD


class aparserAnnotator(ABC):
    @abstractmethod
    def parseTxt(self,s):
        pass
        
### ANNOTATION from SNPEff
class parserAnnotator_ANN(aparserAnnotator):
#ANN             [String]     : Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'     
    def __init__(self, tuplaDictionary):
        descriptionTxt=tuplaDictionary[2].split(':')[1].strip()[1:-1].split('|')
        self.listAnnAnnotation=[]
        for el in descriptionTxt:
            self.listAnnAnnotation.append(el.strip().lower())
        self.POS_MEANING=self.listAnnAnnotation.index('annotation')
        self.POS_ANN_IMPACT=self.listAnnAnnotation.index('annotation_impact')
        self.POS_GENE=self.listAnnAnnotation.index('gene_name')
        self.POS_GENEID=self.listAnnAnnotation.index('gene_id')
        self.POS_FEATURE_TYPE=self.listAnnAnnotation.index('feature_type')
        self.POS_FEATURE_ID=self.listAnnAnnotation.index('feature_id')
        self.POS_TRANSCRIPT=self.listAnnAnnotation.index('transcript_biotype')
        self.POS_AA=self.listAnnAnnotation.index('aa.pos / aa.length')        
        self.POS_HGVSc=self.listAnnAnnotation.index('hgvs.c')        
        self.POS_HGVSp=self.listAnnAnnotation.index('hgvs.p')
        self.POS_RANK=self.listAnnAnnotation.index('rank')
        self.POS_cDNA=self.listAnnAnnotation.index('cdna.pos / cdna.length')
        self.POS_CDS=self.listAnnAnnotation.index('cds.pos / cds.length')       
        self.POS_AA_POSITION=self.listAnnAnnotation.index('aa.pos / aa.length')
    
    def parseTxt(self, s):
    # 0-Allele | 1-Annotation | 2-Annotation_Impact | 3-Gene_Name | 4-Gene_ID | 5-Feature_Type | 6-Feature_ID | 
    # 7-Transcript_BioType | 8-Rank | 9-HGVS.c | 10-HGVS.p | 11-cDNA.pos / cDNA.length | 12-CDS.pos / CDS.length | 
    # 13-AA.pos / AA.length | 14-Distance | 15-ERRORS / WARNINGS / INFO
        affectedGenesSet=[]
        
        # Each string could contains multiple genes/transcript separated by ','
        affectedGenes=s.split(",")
        
        for geneAff in affectedGenes :
            # To each gene/transcript is divided into fields divided by '|'
            fields=geneAff.split('|')
            meaning=fields[self.POS_MEANING]
            impact=fields[self.POS_ANN_IMPACT]
            gene=fields[self.POS_GENE]
            geneID=fields[self.POS_GENEID]
            region=fields[self.POS_FEATURE_TYPE]
            featureID=fields[self.POS_FEATURE_ID]
            transcript=fields[self.POS_TRANSCRIPT]
            AAmodif=fields[self.POS_AA]
            HGVSc=fields[self.POS_HGVSc]
            HGVSp=fields[self.POS_HGVSp]
            
            rank=fields[self.POS_RANK]
            if rank is not None and rank != "" :
                rank=tuple([ int(x) for x in rank.split("/")])
            else :
                rank=None
                
            cdna=fields[self.POS_cDNA]
            if cdna is not None and cdna != "" :
                cdna=tuple([ int(x) for x in cdna.split("/")])
            else :
                cdna=None
            cds=fields[self.POS_CDS]
            if cds is not None and cds != "" :
                cds=tuple([ int(x) for x in cds.split("/")])
            else :
                cds=None
            aa_pos=fields[self.POS_AA_POSITION]
            if aa_pos is not None and aa_pos != "" :
                aa_pos=tuple([ int(x) for x in aa_pos.split("/")])
            else :
                aa_pos=None

            val=typeAnnotation(gene, geneID, meaning, impact, region, featureID, transcript, AAmodif,HGVSc,HGVSp,rank=rank,cDNA=cdna,cds=cds,aa_position=aa_pos)
            
            affectedGenesSet.append(val)
        
        return affectedGenesSet


### ANNOTATION from VEP
class parserAnnotator_CSQ(aparserAnnotator):
#CSQ             [String]     : Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|UNIPROT_ISOFORM|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS 
    def __init__(self, tuplaDictionary):
        descriptionTxt=tuplaDictionary[2].split(':')[1].strip().split('|')
        self.listAnnAnnotation=[]
        for el in descriptionTxt:
            self.listAnnAnnotation.append(el.strip().lower())
        
        self.POS_MEANING=self.listAnnAnnotation.index('consequence')
        self.POS_ANN_IMPACT=self.listAnnAnnotation.index('impact')
        self.POS_GENE=self.listAnnAnnotation.index('symbol')
        self.POS_GENEID=self.listAnnAnnotation.index('gene')
        self.POS_FEATURE_TYPE=self.listAnnAnnotation.index('feature_type')
        self.POS_FEATURE_ID=self.listAnnAnnotation.index('feature')
        self.POS_TRANSCRIPT=self.listAnnAnnotation.index('biotype')
        self.POS_AA=self.listAnnAnnotation.index('amino_acids')
        self.POS_HGVSc=self.listAnnAnnotation.index('hgvsc')
        self.POS_HGVSp=self.listAnnAnnotation.index('hgvsp')
        self.POS_GNOMAD_AF=self.listAnnAnnotation.index('gnomad_af')
            
    def parseTxt(self, s):
        affectedGenesSet=set()
        
        # Each string could contains multiple genes/transcript separated by ','
        affectedGenes=s.split(",")
        
        for geneAff in affectedGenes :
            # To each gene/transcript is divided into fields divided by '|'
            fields=geneAff.split('|')
            meaning=fields[self.POS_MEANING]
            impact=fields[self.POS_ANN_IMPACT]
            gene=fields[self.POS_GENE]
            geneID=fields[self.POS_GENEID]
            region=fields[self.POS_FEATURE_TYPE]
            featureID=fields[self.POS_FEATURE_ID]
            transcript=fields[self.POS_TRANSCRIPT]
            AAmodif=fields[self.POS_AA]
            HGVSc=fields[self.POS_HGVSc]
            HGVSp=fields[self.POS_HGVSp]
            gnomAD_af=None
            try :
                gnomAD_af=float(fields[self.POS_GNOMAD_AF])
            except :
                pass
            val=typeAnnotation(gene, geneID, meaning, impact, region, featureID, transcript, AAmodif,HGVSc,HGVSp,gnomAD_af)
            
            affectedGenesSet.add(val)
        
        return affectedGenesSet

### ANNOTATION from Loss Of Function
class parserAnnotator_LOF(aparserAnnotator):
#LOF             [String]     : Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected' 
    def __init__(self, tuplaDictionary):
        descriptionTxt=tuplaDictionary[2].split(':')[1].strip()[1:-1].split('|')
        self.listAnnAnnotation=[]
        for el in descriptionTxt:
            self.listAnnAnnotation.append(el.strip().lower())
        
        self.POS_GENE=self.listAnnAnnotation.index('gene_name')
        self.POS_GENEID=self.listAnnAnnotation.index('gene_id')
        self.POS_NUM_OF_TRANSCRIPT=self.listAnnAnnotation.index('number_of_transcripts_in_gene')
        self.POS_PERCENT_OF_TRANSCRIPT=self.listAnnAnnotation.index('percent_of_transcripts_affected')
            
    def parseTxt(self, s):
        fields=s[1:-1].split('|')
        gene=fields[self.POS_GENE]
        geneID=fields[self.POS_GENEID]
        try:
            num_transcript=int(fields[self.POS_NUM_OF_TRANSCRIPT])
        except:
            num_transcript=-1
        try:
            percent_transcript=float(fields[self.POS_PERCENT_OF_TRANSCRIPT])
        except:
            percent_transcript=-1
        #print("\t ".join(fields[7:14]))
        return typeAnnotation_LOF(gene, geneID, num_transcript, percent_transcript)

