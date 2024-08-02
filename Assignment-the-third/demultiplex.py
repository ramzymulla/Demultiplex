#!/usr/bin/env python

import bioinfo as bi
import argparse as arp
import math as mth
import gzip
import numpy as np
import matplotlib.pyplot as plt
import itertools as itt

def get_args():
    '''Parses user inputs from command line'''
    parser = arp.ArgumentParser(description="blah blah blah")                   
    parser.add_argument("-q", "--q_cutoff", help="defines Q-score cutoff",         
                        required=True)
    parser.add_argument("-1","--fname1", help="file name for read 1",
                        required=True)
    parser.add_argument("-2","--fname2", help="file name for read 2",
                        required=True)
    parser.add_argument("-3","--fname3", help="file name for read 3",
                        required=True)
    parser.add_argument("-4","--fname4", help="file name for read 4",
                        required=True)
    parser.add_argument("-i","--indices", help="full path to indexes.txt",
                        required=True)
    parser.add_argument("-p", "--fpath", help="path to file directory",
                        required=True)
    parser.add_argument("-o", "--output", help="path to output directory",
                        required=True)
    return parser.parse_args()

def revcomp(seq: str) -> str:
    rc = ''
    for c in seq[::-1]:
        rc+=bi.DNA_COMPS[c]
    return rc

def get_record(file) -> tuple:
    header=file.readline().strip()
    if header=='':
        return 0,0,0
    seq = file.readline().strip()
    file.readline()
    qscore=file.readline().strip()
    return header,seq,qscore

def get_indices(file:str)->list:
    barcodes = []
    with open(file,'r') as f:
        line = f.readline().strip().split()
        while line:
            barcodes.append(line[-1])
            line = f.readline().strip().split()
    return barcodes

def check_qscore(qscores: str,cutoff:int)->bool:
    mean = np.mean([bi.qual_score(i) for i in qscores])
    return bool(mean > cutoff)

args = get_args()
q_cutoff = args.q_cutoff
path = args.fpath
outpath = args.output
barcodes = get_indices(args.indices)
headers,seqs,qscores={},{},{}
matches,unknowns,ihops = 0,0,0
index_pairs={}
known_pairs=[]

for i in barcodes:
    for j in barcodes:
        known_pairs.append(f"{i}-{j}")
        index_pairs[f"{i}-{j}"] = 0


r1_outs = {f"{bar}-{bar}":gzip.open(f"{outpath}read1_{bar}-{bar}.fastq",'wt') for bar in barcodes}
r2_outs = {f"{bar}-{bar}":gzip.open(f"{outpath}read2_{bar}-{bar}.fastq",'wt') for bar in barcodes}
r1_outs["unknowns"] = gzip.open(f"{outpath}read1_unknowns.fastq",'wt')
r2_outs["unknowns"] = gzip.open(f"{outpath}read2_unknowns.fastq",'wt')
r1_outs["ihops"] = gzip.open(f"{outpath}read1_ihops.fastq",'wt')
r2_outs["ihops"] = gzip.open(f"{outpath}read2_ihops.fastq",'wt')



with gzip.open(path+args.fname1,'rt') as R1,\
gzip.open(path+args.fname2,'rt') as R2,\
gzip.open(path+args.fname3,'rt') as R3,\
gzip.open(path+args.fname4,'rt') as R4:
    
    handles = [R1,R2,R3,R4]
    
    while True:
        hr1,sr1,qr1=get_record(R1)
        if hr1==0: break
        hi1,si1,qi1=get_record(R2)
        hi2,si2_rc,qi2_r=get_record(R3)
        hr2,sr2,qr2=get_record(R4)

        si2 = revcomp(si2_rc)
        qi2 = qi2_r[::-1]

        bar = f"{si1}-{si2}"

        h1 = f"{hr1} {bar}"
        h4 = f"{hr2} {bar}"

        record1 = f"{hr1}\n{sr1}\n{qr1}"
        record2 = f"{hr2}\n{sr2}\n{qr2}"

        if (bar not in known_pairs) or\
            (not check_qscore(qi1,q_cutoff)) or\ # type: ignore
            (not check_qscore(qi2,q_cutoff)): # type: ignore
            unknowns+=1
            if bar not in index_pairs:
                index_pairs[bar]=1
            else:
                index_pairs[bar]+=1
            if unknowns ==1: 
                r1_outs['unknowns'].write(record1)
                r2_outs['unknowns'].write(record2)
            else:
                r1_outs['unknowns'].write("\n"+record1)
                r2_outs['unknowns'].write("\n"+record2)
            
        elif si1==si2:
            matches+=1
            index_pairs[bar]+=1
            if matches==1:
                r1_outs[bar].write(record1)
                r2_outs[bar].write(record2)
            else:
                r1_outs['matches'].write("\n"+record1)
                r2_outs['matches'].write("\n"+record2)
        else:
            ihops+=1
            index_pairs[bar]+=1
            if ihops==1:
                r1_outs['ihops'].write(record1)
                r2_outs['ihops'].write(record2)
            else:
                r1_outs['ihops'].write("\n"+record1)
                r2_outs['ihops'].write("\n"+record2)

        # header1 = R1.readline().strip() 
        # if header1 == "": break
        # header2 = R4.readline().strip()
        # read1 = R1.readline().strip()
        # read2 = R4.readline().strip()
        # R2.readline()
        # R3.readline()
        # bar1 = R2.readline().strip()
        # bar2_rc = R3.readline().strip()
        # bar2 = revcomp(bar2_rc) # type: ignore
        # fullbar = f"{bar1}-{bar2}"
        # header1 = f"{header1} {fullbar}"
        # header2 = f"{header2} {fullbar}"
        # for f in handles:
        #     f.readline()
        # qual1 = R2.readline().strip()
        # qual2 = R3.readline().strip()
        # R2.readline()
        # R3.readline()
        # record1 = f"{header1}\n{read1}\n+\n{qual1}"
        # record2 = f"{header2}\n{read2}\n+\n{qual2}"

        # if (fullbar not in known_pairs) or\
        #     (not check_qscore(qual1,q_cutoff)) or\ # type: ignore
        #     (not check_qscore(qual2,q_cutoff)): # type: ignore
        #     unknowns+=1
        #     if fullbar not in index_pairs:
        #         index_pairs[fullbar]=1
        #     else:
        #         index_pairs[fullbar]+=1
        #     if unknowns ==1: 
        #         r1_outs['unknowns'].write(record1)
        #         r2_outs['unknowns'].write(record2)
        #     else:
        #         r1_outs['unknowns'].write("\n"+record1)
        #         r2_outs['unknowns'].write("\n"+record2)
            
        # elif bar1==bar2:
        #     matches+=1
        #     index_pairs[fullbar]+=1
        #     if matches==1:
        #         r1_outs[fullbar].write(record1)
        #         r2_outs[fullbar].write(record2)
        #     else:
        #         r1_outs['matches'].write("\n"+record1)
        #         r2_outs['matches'].write("\n"+record2)
        # else:
        #     ihops+=1
        #     index_pairs[fullbar]+=1
        #     if ihops==1:
        #         r1_outs['ihops'].write(record1)
        #         r2_outs['ihops'].write(record2)
        #     else:
        #         r1_outs['ihops'].write("\n"+record1)
        #         r2_outs['ihops'].write("\n"+record2)

for f in r1_outs: r1_outs[f].close()
for f in r2_outs: r2_outs[f].close()