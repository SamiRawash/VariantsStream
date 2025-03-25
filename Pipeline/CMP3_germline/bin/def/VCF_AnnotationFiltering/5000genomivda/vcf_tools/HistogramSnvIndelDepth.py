#! /usr/bin/env python
import argparse, os, numpy as np #, math
#import numba as nb

ACTIVATE_NUMBA=True # False # 
COLUMN_FORMAT=8
COLUMN_GENOTYPE=9
        
n_depth=0
n_af=0

depth_min=0
depth_max=200
depth_step=5

af_min=0
af_max=1
af_step=0.01


#@nb.jit(nopython=ACTIVATE_NUMBA)
#@guvectorize([(float64[:], float64[:])], '()->()')
def extractAD_and_AF(ad):
    splitted=ad.split(",")
    alt=int(splitted[1])
    ref=int(splitted[0])
    cov=ref+alt
    #print(f"{alt}\t{ref}")
    return cov, alt, alt/cov if cov>0 else -1

#@nb.jit(nopython=ACTIVATE_NUMBA)
def parseRow(mydata):
    
    depth=np.zeros(n_depth+1, dtype=int) # [ 0 ] * (n_depth +1 ) #
    af=np.zeros(n_af+1, dtype=int) # [ 0 ] * (n_af +1 ) # 
    alt=np.zeros(n_depth+1, dtype=int) # [ 0 ] * (n_depth +1 ) #
    
    for row in mydata :
        f=row[0].split(":")
        g=row[1].split(":")
        AD=g[f.index("AD")]
        
        c,d,a,=extractAD_and_AF(str(AD))
        if c>=depth_min and c<depth_max :
            depth[int(np.floor_divide(c-depth_min,depth_step))]+=1
        else :
            depth[n_depth]+=1
            
        if d>=depth_min and d<depth_max :
            alt[int(np.floor_divide(d-depth_min,depth_step))]+=1
        else :
            alt[n_depth]+=1
        
        if a>=af_min and a<af_max :
            af[int(np.floor_divide(a-af_min,af_step))]+=1
        else :
            af[n_af]+=1
    
    return depth,alt,af

"""
#@nb.jit(nopython=ACTIVATE_NUMBA)
def parseline(line):
    splitted=line.strip().split("\t")
    f=splitted[COLUMN_FORMAT].split(":")
    g=splitted[COLUMN_GENOTYPE].split(":")
    AD=g[f.index("AD")]
    AD=[ int(x) for x in AD.split(",") ]
    return AD[1],float(AD[1])/float((AD[0]+AD[1]))

#@nb.jit(nopython=ACTIVATE_NUMBA)
def readfile(file):
    line=file.readline()
    while line != "" :
        if line[0] != "#" :
            
            d,a=parseline(line)
            
            if d>=depth_min and d<depth_max :
                depth[int(np.floor_divide(d-depth_min,depth_step))]+=1
            else :
                depth[n_depth]+=1
                
            if a>=af_min and a<af_max :
                af[int(np.floor_divide(a-af_min,af_step))]+=1
            else :
                af[n_af]+=1
        line=file.readline()        
    return 0


class HistogramDepth():
    
    @nb.jit(nopython=ACTIVATE_NUMBA)
    def parseline(self,line):
        splitted=line.strip().split("\t")
        f=splitted[self.COLUMN_FORMAT].split(":")
        g=splitted[self.COLUMN_GENOTYPE].split(":")
        AD=g[f.index("AD")]
        AD=[ int(x) for x in AD.split(",") ]
        return AD[1],float(AD[1])/float((AD[0]+AD[1]))
    
    @nb.jit(nopython=ACTIVATE_NUMBA)
    def readfile(self,file):
        line=file.readline()
        while line != "" :
            if line[0] != "#" :
                
                d,a=self.parseline(line)
                
                if d>=self.depth_min and d<self.depth_max :
                    self.depth[int(np.floor_divide(d-self.depth_min,self.depth_step))]+=1
                else :
                    self.depth[self.n_depth]+=1
                    
                if a>=self.af_min and a<self.af_max :
                    self.af[int(np.floor_divide(a-self.af_min,self.af_step))]+=1
                else :
                    self.af[self.n_af]+=1
            line=file.readline()        
        return 0
            
    def __init__(self, vcf, depth_min, depth_max, depth_step, af_min, af_max, af_step, genotype_column=9):
        self.COLUMN_FORMAT=8
        self.n_depth=int((depth_max-depth_min)/depth_step)
        self.n_af=int((depth_max-depth_min)/depth_step) 
        self.COLUMN_GENOTYPE=genotype_column
        
        self.depth_min=depth_min
        self.depth_max=depth_max
        self.depth_step=depth_step
        self.depth=np.zeros(self.n_depth+1, dtype=int)
        
        self.af_min=af_min
        self.af_max=af_max
        self.af_step=af_step
        self.af=np.zeros(self.n_af+1, dtype=int)
        with open(vcf, "r") as file :
            self.readfile(file)
"""
    
################################################################
#### MAIN
################################################################ 
if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="Produce Depth Histograms from VCF: ")
    parser.add_argument('-vcf','--vcf',action='store',type=str,help="VCF source data file.", required=True, default=None)
    parser.add_argument('-depth_min','--depth_min',action='store',type=int,help="Minimum depth of coverage in the histogram.", required=False, default=0)
    parser.add_argument('-depth_max','--depth_max',action='store',type=int,help="Maximum depth of coverage in the histogram.", required=False, default=500)
    parser.add_argument('-depth_step','--depth_step',action='store',type=int,help="Maximum depth of coverage in the histogram.", required=False, default=1)
    parser.add_argument('-af_min','--af_min',action='store',type=float,help="Minimum allele frequency of coverage in the histogram.", required=False, default=0.0)
    parser.add_argument('-af_max','--af_max',action='store',type=float,help="Maximum allele frequency of coverage in the histogram.", required=False, default=1.0)
    parser.add_argument('-af_step','--af_step',action='store',type=float,help="Maximum allele frequency of coverage in the histogram.", required=False, default=0.01)
    parser.add_argument('-c','--genotype_column',action='store',type=float,help="Column.", required=False, default=9)
    parser.add_argument('-od','--outdp',action='store',type=str,help="Histogram data of the coverage depth.", required=False, default="OutHistCoverageDepth.tsv")
    parser.add_argument('-ol','--outal',action='store',type=str,help="Histogram data of the alteration allel depth.", required=False, default="OutHistDepthAllel.tsv")
    parser.add_argument('-oa','--outaf',action='store',type=str,help="Histogram data of the alteration allel frequency.", required=False, default="OutHistAllelFreq.tsv")
    arg=parser.parse_args()
    
    COLUMN_GENOTYPE=arg.genotype_column
    
    depth_min=arg.depth_min
    depth_max=arg.depth_max
    depth_step=arg.depth_step
    
    af_min=arg.af_min
    af_max=arg.af_max
    af_step=arg.af_step
            
    n_depth=int((depth_max-depth_min)/depth_step)
    n_af=int((af_max-af_min)/af_step)
    """
    depth=np.zeros(n_depth+1, dtype=int)
    af=np.zeros(n_af+1, dtype=int)
    """
    file_list=arg.vcf.split(",")
    depth=np.zeros((len(file_list)+2,n_depth+1), dtype=float) # [ 0 ] * (n_depth +1 ) #
    alt=np.zeros((len(file_list)+2,n_depth+1), dtype=float)
    af=np.zeros((len(file_list)+2,n_af+1), dtype=float) 
    for ndx,file in  enumerate(file_list) :
        data=np.loadtxt(file, usecols=[COLUMN_FORMAT,COLUMN_GENOTYPE], dtype=str)
        depth_cc,alt_cc,af_cc=parseRow(data)
        depth[ndx+2,:]=depth_cc
        alt[ndx+2,:]=alt_cc
        af[ndx+2,:]=af_cc
    
    header="#Min\tMax\t"+"\t".join([os.path.basename(x) for x in file_list])

    min_column=np.arange(depth_min,depth_max+depth_step,depth_step, dtype=float)
    depth[0,:-1]=min_column[:-1]
    depth[1,:-1]=min_column[1:]
    depth[0,-1]=np.NaN
    depth[1,-1]=np.NaN
    
    column_format=['%-.0f','%-.0f']
    column_format.extend(['%.0f']*len(file_list))
    np.savetxt(arg.outdp,depth.T, fmt=tuple(column_format), delimiter='\t', newline='\n', header=header)
    
    alt[0,:-1]=min_column[:-1]
    alt[1,:-1]=min_column[1:]
    alt[0,-1]=np.NaN
    alt[1,-1]=np.NaN
    np.savetxt(arg.outal,alt.T, fmt=tuple(column_format), delimiter='\t', newline='\n', header=header)


    min_column=np.arange(af_min,af_max+af_step,af_step)
    af[0,:-1]=min_column[:-1]
    af[1,:-1]=min_column[1:]
    af[0,-1]=np.NaN
    af[1,-1]=np.NaN
    
    column_format=['%-.3f','%-.3f']
    column_format.extend(['%.0f']*len(file_list))
    np.savetxt(arg.outaf,af.T, fmt=tuple(column_format), delimiter='\t', newline='\n', header=header)
    

    """
    with open(arg.vcf, "r") as file :
        readfile(file)
    
    #historgam=HistogramDepth(,)
    #print(f"Depth Avg: {np.average(depth)} Std: {np.std(depth)}")
    with open(arg.outdp, "w") as outf :
        outf.write("#{:7}\t{:8}\t{:8}\n".format("Min","Max","Count"))
        for ndx,el in enumerate(depth) :
            if ndx == n_depth :
                outf.write("{:8}\t{:8}\t{:<d}\n".format("Nan","Nan",el))
            else :
                outf.write("{:<8d}\t{:<8d}\t{:<d}\n".format(depth_min+ndx*depth_step, depth_min+(ndx+1)*depth_step, el))
    
    with open(arg.outaf, "w") as outf :
        outf.write("#{:7}\t{:8}\t{:8}\n".format("Min","Max","Count"))
        for ndx,el in enumerate(af) :
            if ndx == n_af :
                outf.write("{:8}\t{:8}\t{:<d}\n".format("Nan","Nan",el))
            else :
                outf.write("{:<8.3f}\t{:<8.3f}\t{:<d}\n".format(af_min+ndx*af_step, af_min+(ndx+1)*af_step, el))
    """
    