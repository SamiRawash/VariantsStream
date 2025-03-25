#! /usr/bin/env python3

"""
class score_predictor():
    
    def __init__(self):
        self.value_type=None
        pass
    
    def pathogenicity(self, value):
        return False
    
class score_predictor_FATHMM_score(score_predictor):
    
    def __init__(self, cutoff_value=-1.5):
        self.value_type=float
        self.cutoff_value=cutoff_value
    
    def pathogenicity(self, value):
        if isinstance(value, self.value_type):
            if value<=self.cutoff_value:
                return True
            else:
                return False
    
"""
PATHOGENICITY_SCORE_HEADER=tuple(['Pathogenic', 'Benign', 'Uncertain,Unknown'])

def score_predictor_FATHMM_score(value):
    cutoff_value=-1.5
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value<=cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError:
            return (0,0,1)
    else:
        return (0,0,1)

# https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
def score_predictor_CLINVAR_score(value):
    set_benign=set(["benign","benign/likely_benign","likely_benign"])
    set_Uncertain=set(["uncertain_significance"])
    if value is not None and isinstance(value, str):
        new_value=value.lower()
        if new_value in set_benign:
            return (0,1,0)
        elif new_value in set_Uncertain:
            return (0,0,1)
        else:
            return (1,0,0)
    else:
        return (0,0,1) 
    
def score_predictor_SIFT_score(value):
    cutoff_value=0.05
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value<=cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError:
            return (0,0,1)
    else:
        return (0,0,1)
    
def score_predictor_POLYPHEN2_HDIV_score(value):
    cutoff_value=0.446
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value>cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError:
            return (0,0,1)
    else:
        return (0,0,1)

def score_predictor_POLYPHEN2_HVAR_score(value):
    cutoff_value=0.446
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value>cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError:
            return (0,0,1)
    else:
        return (0,0,1)

def score_predictor_CADD_score(value):
    cutoff_value=30
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value>=cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError:
            return (0,0,1)
    else:
        return (0,0,1)

def score_predictor_DANN_score(value):
    cutoff_value=0.90
    if value is not None :
        try: 
            if not isinstance(value, float):
                value=float(value)
            if value>=cutoff_value:
                return (1,0,0)
            else:
                return  (0,1,0)
        except ValueError :
            return (0,0,1)
    else:
        return (0,0,1)
    
scorepredictionDict={}
scorepredictionDict['FATHMM_score']=score_predictor_FATHMM_score #score_predictor_FATHMM_score()
scorepredictionDict['CLNSIG']=score_predictor_CLINVAR_score
scorepredictionDict['SIFT_score']=score_predictor_SIFT_score
scorepredictionDict['Polyphen2_HDIV_score']=score_predictor_POLYPHEN2_HDIV_score
scorepredictionDict['Polyphen2_HVAR_score']=score_predictor_POLYPHEN2_HVAR_score
scorepredictionDict['CADD_phred']=score_predictor_CADD_score
scorepredictionDict['DANN_score']=score_predictor_DANN_score

if __name__ == '__main__' :
    info={'FATHMM_score':'-2', 'CLNSIG':'Benign', 'SIFT_score':'0.03', 'Polyphen2_HDIV_score':'0.222', 'Polyphen2_HVAR_score':'0.698', 'CADD_phred':'31'}
    
    set_score_pred=scorepredictionDict.keys()
    counter_pathogenic=0
    counter_benign=0
    counter_uncertain=0
    counter_predictor=0
    
    for key_score_pred in set_score_pred :
        if key_score_pred in info :
            counter_predictor+=1
            sp=scorepredictionDict[key_score_pred]
            res=sp(info[key_score_pred])
            counter_pathogenic+=res[0]
            counter_benign+=res[1]
            counter_uncertain+=res[2]
            if res[0]==1:
                print("pathogenic")
            else:
                print("non pathogenic")
"""
if 'FATHMM_score' in info and 'FATHMM_score' in scorepredictionDict:
    sp=scorepredictionDict['FATHMM_score']
    res=sp(info['FATHMM_score'])
    if res:
        counter+=1
        print("pathogenic")
    else:
        print("non pathogenic")
"""
