#!/usr/bin/env python

import bioinfo as bi
import argparse as arp
import gzip
import numpy as np
import matplotlib.pyplot as plt

DNA_TAB = str.maketrans(bi.DNA_COMPS)
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
    '''
    Converts input DNA sequence to its reverse complement

    Args:
        seq (str): DNA sequence

    Returns:
        str: reverse complement of seq
    '''
    return seq[::-1].translate(DNA_TAB)

def get_record(file) -> tuple:
    '''
    FASTQ file handle

    Args:
        file (io.TextWrapper): open FASTQ file handle (at line x such that
                                x+1 is a header line)

    Returns:
        tuple: header, sequence, qscore 
        **file line pointer now at x+4
    '''
    header=file.readline().strip()
    if header=='':
        return 0,0,0
    seq = file.readline().strip()
    file.readline()
    qscore=file.readline().strip()
    return header,seq,qscore

def get_indices(file:str)->list:
    '''
    Extracts indices/barcodes from tab separated indexes.txt

    Args:
        file (str): tab separated text file containing a list of 
                    known indices/barcodes (indices must be in last column)

    Returns:
        list: list of known indices/barcodes
    '''
    barcodes = []
    with open(file,'r') as f:
        line = f.readline().strip().split()
        line = f.readline().strip().split()
        while line:
            barcodes.append(line[-1])
            line = f.readline().strip().split()
    return barcodes

def check_qscore(qscores: str,cutoff:int,enc=33)->bool:
    '''
    Checks if the qscore is above the cutoff value

    Args:
        qscores (str): string of qscores
        cutoff (int): desired cutoff value
        phred (int): phred encoding (default is 33)

    Returns:
        bool: True if avg qscore >= cutoff
    '''
    num = 0
    for i in qscores:
        num += bi.convert_phred(i,offset=enc)
    return num/len(qscores) >= cutoff

args = get_args()
q_cutoff = int(args.q_cutoff)
path = args.fpath
outpath = args.output
barcodes = get_indices(args.indices)

matches,unknowns,ihops,lowq = 0,0,0,0
ind_revcomps={revcomp(bar):bar for bar in barcodes}
all_pairs={}
match_pairs={}
ihopped_pairs={}

for i in barcodes:
    for j in barcodes:
        all_pairs[f"{i}-{j}"] = 0
        if i==j:
            match_pairs[f"{i}-{j}"] = 0

r1_outs = {bar:open(f"{outpath}read1_{bar}.fastq",'w') for bar in match_pairs}
r2_outs = {bar:open(f"{outpath}read2_{bar}.fastq",'w') for bar in match_pairs}
r1_outs["unknowns"] = open(f"{outpath}read1_unknowns.fastq",'w')
r2_outs["unknowns"] = open(f"{outpath}read2_unknowns.fastq",'w')
r1_outs["ihops"] = open(f"{outpath}read1_ihops.fastq",'w')
r2_outs["ihops"] = open(f"{outpath}read2_ihops.fastq",'w')

# with open(path+args.fname1,'r') as R1,\
# open(path+args.fname2,'r') as R2,\
# open(path+args.fname3,'r') as R3,\
# open(path+args.fname4,'r') as R4:

with gzip.open(path+args.fname1,'rt') as R1,\
gzip.open(path+args.fname2,'rt') as R2,\
gzip.open(path+args.fname3,'rt') as R3,\
gzip.open(path+args.fname4,'rt') as R4:
    while True:
        hr1,sr1,qr1=get_record(R1)
        if hr1==0: break
        hi1,si1,qi1=get_record(R2)
        hi2,si2_rc,qi2=get_record(R3)
        hr2,sr2,qr2=get_record(R4)
        if si2_rc in ind_revcomps:
            si2 = ind_revcomps[si2_rc]
        else:
            si2 = revcomp(si2_rc)
        # qi2 = qi2[::-1]
        bar = f"{si1}-{si2}"
        # print(bar)

        hn1 = f"{hr1} {bar}"
        hn2 = f"{hr2} {bar}"

        record1 = f"{hn1}\n{sr1}\n+\n{qr1}"
        record2 = f"{hn2}\n{sr2}\n+\n{qr2}"
        
        if bar in all_pairs: 
            # if not check_qscore(qi1+qi2,q_cutoff):
            if (not check_qscore(qi1,q_cutoff)) or (not check_qscore(qi2,q_cutoff)):
                lowq += 1
                unknowns += 1
                if unknowns ==1: 
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs['unknowns'].write(record1)
                    r2_outs['unknowns'].write(record2)
                else:
                    r1_outs['unknowns'].write("\n"+record1)
                    r2_outs['unknowns'].write("\n"+record2)
                    
            elif bar in match_pairs:
                matches+=1
                match_pairs[bar]+=1
                if match_pairs[bar]==1:
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs[bar].write(record1)
                    r2_outs[bar].write(record2)
                else:
                    r1_outs[bar].write("\n"+record1)
                    r2_outs[bar].write("\n"+record2)

            else:
                ihops+=1
                if bar in ihopped_pairs:
                    ihopped_pairs[bar]+=1
                else:
                    ihopped_pairs[bar]=1
                if ihops==1:
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs['ihops'].write(record1)
                    r2_outs['ihops'].write(record2)
                else:
                    r1_outs['ihops'].write("\n"+record1)
                    r2_outs['ihops'].write("\n"+record2)
        else:
            unknowns+=1
            if unknowns ==1: 
                # print(record1+"\n"+record2+"\n\n")
                r1_outs['unknowns'].write(record1)
                r2_outs['unknowns'].write(record2)
            else:
                r1_outs['unknowns'].write("\n"+record1)
                r2_outs['unknowns'].write("\n"+record2)

for f in r1_outs: r1_outs[f].close()
for f in r2_outs: r2_outs[f].close()
# print(f"Number of Index-Hopped Reads: {ihops}")
# print(f"Number of Unknown Reads: {unknowns}")
# for i in match_pairs:
#     print(f"Number of Reads for barcode {i}: {match_pairs[i]}")
# for i in index_pairs:
#     if i not in match_pairs and index_pairs[i] != 0:
#         print(f"Number of Reads with hopped barcode {i}: {index_pairs[i]}")
total_reads=matches+ihops+unknowns
f = open(f"{outpath}/stats.txt",'w')
f.write(f"Barcode_1\tBarcode_2\tNum_Reads\tPercentage")
for bar in match_pairs:
    bar1,bar2 = bar.split('-')
    val=match_pairs[bar]
    percentage=np.round(100*(val/total_reads),2)
    f.write(f"\n{bar1}\t{bar2}\t{val}\t{percentage}")
for bar in ihopped_pairs:
    bar1,bar2 = bar.split('-')
    val=ihopped_pairs[bar]
    percentage=np.round(100*(val/total_reads),2)
    f.write(f"\n{bar1}\t{bar2}\t{val}\t{percentage}")
f.close()
print(f"Total Matched Reads: {matches} ({np.round(100*(matches/total_reads),2)}%)")
print(f"Total Index-Hopped Reads: {ihops} ({np.round(100*(ihops/total_reads),2)}%)")
print(f"Total Unknown Reads: {unknowns} ({np.round(100*(unknowns/total_reads),2)}%)")
print(f"Total Reads Below Cutoff: {lowq} ({np.round(100*(lowq/total_reads),2)}%)")
