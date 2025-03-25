import argparse

# define the input and output file paths
parser = argparse.ArgumentParser(description="vcf SV")
parser.add_argument("--input",action="store",type=str,help="input file")
parser.add_argument("--sample_id",action="store",type=str,help="sample id of consensus")
arg=parser.parse_args()

with open(arg.input, 'r') as vcf, \
    open(arg.sample_id+"_DEL.vcf","w") as del_vcf, \
    open(arg.sample_id+"_DUP.vcf","w") as dup_vcf, \
    open(arg.sample_id+"_INV.vcf","w") as inv_vcf, \
    open(arg.sample_id+"_TRA.vcf","w") as tra_vcf:
    
    for line in vcf.readlines():
        if line[0]=="#":
            del_vcf.write(line)
            dup_vcf.write(line)
            inv_vcf.write(line)
            tra_vcf.write(line)
        elif "SVTYPE=DEL" in line:
            del_vcf.write(line)
        elif "SVTYPE=DUP" in line:
            dup_vcf.write(line)
        elif "SVTYPE=INV" in line:
            inv_vcf.write(line)
        elif "SVTYPE=TRA" in line:
            tra_vcf.write(line)

del_vcf.close()
dup_vcf.close()
inv_vcf.close()
tra_vcf.close()