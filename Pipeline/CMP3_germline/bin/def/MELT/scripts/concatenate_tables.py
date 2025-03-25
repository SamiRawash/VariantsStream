#! /usr/bin/env python3
import pandas, argparse, os

if(__name__=='__main__'):
    parser = argparse.ArgumentParser(description="This script takes in input a list of tables (nice to give the complete path) and concatenate in one unique table.")
    ### INPUT - File
    parser.add_argument('-i','--input',action='store',type=str,help="List of files you want to concatenate, comma separated.", required=True, default=None)
    parser.add_argument('-o', '--output', action='store',help="Output file.", type=str, required=True, default=None)
    arg=parser.parse_args()


    list_files=arg.input
    out=arg.output

    file=list_files.split(",")
    thelist=[]
    for filename in file:
        if os.path.exists(filename) == True:
            df = pandas.read_csv(filename, sep="\t")
            thelist.append(df)
    joined=pandas.concat(thelist, axis=0, ignore_index=True)
    joined.to_csv(out, sep ="\t", index=False)
