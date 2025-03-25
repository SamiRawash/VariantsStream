'''
Created on 9 Mar 2023

@author: flanduzzi
'''
def toFloat(vect):
    x=[None]*len(vect)
    for ndx,v in enumerate(vect) :
        try: 
            x[ndx]=float(v)
        except: pass
    return x

def generate_new_filter(convert,strong_ben,mod_ben,supp_ben,supp_pat,mod_pat,strong_pat):
    
    def my_criteria(x):
        if x is None or x==".":
            return 0
        try :
            x=float(x)*convert
            if strong_ben is not None and x<=strong_ben:
                return -3
            elif mod_ben is not None and x<=mod_ben:
                return -2
            elif supp_ben is not None and x<=supp_ben:
                return -1
            elif supp_pat is not None and x>=supp_ben:
                return 1
            elif mod_pat is not None and x>=mod_pat:
                return 2
            elif strong_pat is not None and x>=strong_pat:
                return 3
        except : 
            pass
        return 0
        
    return lambda x: my_criteria(x)
    

def loadVarsomeCriteriaTable(file_path):
    dictCriteria={}
    with open(file_path) as f :
        lines=f.readlines()
        for line in lines :
            if len(line)>0 and line[0]!='#' :
                splitted=line.split('\t')
                identifier=splitted[1].strip()
                if identifier is not None and identifier!="" :
                    name=splitted[0].strip()
                    conv_factor=int(splitted[2])
                    calibration=int(splitted[3])
                    vals=toFloat(splitted[3:])
                    dictCriteria[identifier]=(name,conv_factor,calibration,vals, generate_new_filter(conv_factor,vals[4],vals[5],vals[6],vals[8],vals[9],vals[10]))
    return dictCriteria        


