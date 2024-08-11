#!/usr/bin/env python

import bioinfo as bi
import argparse as arp
import gzip
import numpy as np
import matplotlib.pyplot as plt

DNA_TAB = str.maketrans(bi.DNA_COMPS)           # translation table for DNA complements
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
    parser.add_argument("-p", "--fpath", help="path to input FASTQ directory",
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
    return seq[::-1].translate(DNA_TAB)      # returns translation of reversed sequence

def get_record(file) -> tuple:
    '''
    FASTQ file handle

    Args:
        file (io.TextWrapper): FASTQ file handle 
        (NOTE: must be at line x such that x+1 is a header line)

    Returns:
        tuple: header, sequence, qscore 
        **file line pointer now at x+4
    '''
    header=file.readline().strip()          # extract head line
    if header=='':                          # returns 0's if EOF
        return 0,0,0
    seq = file.readline().strip()           # extract sequence line
    file.readline()                         # skip "+" line
    qscore=file.readline().strip()          # extract qscore line

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
    barcodes = []                                   # initialize list
    with open(file,'r') as f:     
        line = f.readline().strip().split()         # skip header line
        line = f.readline().strip().split()         # load in first row
        while line:
            barcodes.append(line[-1])               # append last column to barodes
            line = f.readline().strip().split()     # iterate to next row until EOF
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
    num = 0                                     # initialize accumulator
    for i in qscores:
        num += bi.convert_phred(i,offset=enc)   # convert each character and add to num
    return num/len(qscores) >= cutoff

args = get_args()                               # get input args
q_cutoff = int(args.q_cutoff)                   # convert q_cutoff to int
path = args.fpath                               # assign args to regular variables
outpath = args.output
barcodes = get_indices(args.indices)            # extract barcodes

matches,unknowns,ihops,lowq = 0,0,0,0           # initialize accumulators

ind_revcomps={revcomp(bar):bar for bar in barcodes}     # make revcomps dict
all_pairs={}                                    # dict for all possible pairs
match_pairs={}                                  # dict for only matched pairs
ihopped_pairs={}                                # dict for index-hopped pairs

for i in barcodes:
    for j in barcodes:
        all_pairs[f"{i}-{j}"] = 0               # initialize each entree to 0
        if i==j:
            match_pairs[f"{i}-{j}"] = 0
            
### open write files into r1_ and r2_ outs dicts
r1_outs = {bar:open(f"{outpath}read1_{bar}.fastq",'w') for bar in match_pairs}
r2_outs = {bar:open(f"{outpath}read2_{bar}.fastq",'w') for bar in match_pairs}
r1_outs["unknowns"] = open(f"{outpath}read1_unknowns.fastq",'w')
r2_outs["unknowns"] = open(f"{outpath}read2_unknowns.fastq",'w')
r1_outs["ihops"] = open(f"{outpath}read1_ihops.fastq",'w')
r2_outs["ihops"] = open(f"{outpath}read2_ihops.fastq",'w')

# with open(path+args.fname1,'r') as R1,\ 
# open(path+args.fname2,'r') as R2,\
# open(path+args.fname3,'r') as R3,\
# open(path+args.fname4,'r') as R4:                     # normal input files (unit tests)


with gzip.open(path+args.fname1,'rt') as R1,\
gzip.open(path+args.fname2,'rt') as R2,\
gzip.open(path+args.fname3,'rt') as R3,\
gzip.open(path+args.fname4,'rt') as R4:                 # zipped input files
    
    while True:
        hr1,sr1,qr1=get_record(R1)          # get read 1 record
        if hr1==0: break                    # break loop if EOF
        hi1,si1,qi1=get_record(R2)          # get index 1 record
        hi2,si2_rc,qi2=get_record(R3)       # get index 2 record
        hr2,sr2,qr2=get_record(R4)          # get read 2 record

        if si2_rc in ind_revcomps:          # reverse index 2 sequence
            si2 = ind_revcomps[si2_rc]      # use dict if known index
        else:
            si2 = revcomp(si2_rc)           # else use revcomp()
        # qi2 = qi2[::-1]
        bar = f"{si1}-{si2}"                # make barcode
        # print(bar)

        hn1 = f"{hr1} {bar}"                # make new header lines
        hn2 = f"{hr2} {bar}"

        record1 = f"{hn1}\n{sr1}\n+\n{qr1}"        # make new records
        record2 = f"{hn2}\n{sr2}\n+\n{qr2}"
        
        if bar in all_pairs:                       # start with matched and i-hopped
            if (not check_qscore(qi1,q_cutoff)) or \
                (not check_qscore(qi2,q_cutoff)):               # check index q-scores
                lowq += 1                                       # increment lowq and unknowns accumulators
                unknowns += 1
                if unknowns ==1:                                # check if first record
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs['unknowns'].write(record1)          # write to unknowns
                    r2_outs['unknowns'].write(record2)
                else:                                          
                    r1_outs['unknowns'].write("\n"+record1)
                    r2_outs['unknowns'].write("\n"+record2)              
            elif bar in match_pairs:                            # check if matched pairs
                matches+=1                                      # increment accumulators
                match_pairs[bar]+=1
                if match_pairs[bar]==1:                         # check if first record
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs[bar].write(record1)                 # write to sample file
                    r2_outs[bar].write(record2)
                else:
                    r1_outs[bar].write("\n"+record1)
                    r2_outs[bar].write("\n"+record2)
            else:                                               # not matched => i-hopped
                ihops+=1                                        # increment accumulators
                if bar in ihopped_pairs:                        # track occurance of each i-hopped pair
                    ihopped_pairs[bar]+=1
                else:
                    ihopped_pairs[bar]=1                        # initialize to 0 if first time seeing this pair
                
                if ihops==1:                                    # write to i-hops file
                    # print(record1+"\n"+record2+"\n\n")
                    r1_outs['ihops'].write(record1)
                    r2_outs['ihops'].write(record2)
                else:
                    r1_outs['ihops'].write("\n"+record1)
                    r2_outs['ihops'].write("\n"+record2)

        else:                      # all remaining must one or more "N" or SNP
            unknowns+=1                                         # increment accumulator
            if unknowns ==1:                                    # write to unknowns file
                # print(record1+"\n"+record2+"\n\n")
                r1_outs['unknowns'].write(record1)
                r2_outs['unknowns'].write(record2)
            else:
                r1_outs['unknowns'].write("\n"+record1)
                r2_outs['unknowns'].write("\n"+record2)

for f in r1_outs: r1_outs[f].close()                            # close output FASTQs
for f in r2_outs: r2_outs[f].close()

total_reads=matches+ihops+unknowns                              # calculate total read count
f = open(f"{outpath}/stats.txt",'w')                            # open stats file (tab delim)
f.write(f"Barcode_1\tBarcode_2\tNum_Reads\tPercentage")         # write header line
for bar in match_pairs:                                         # write stats for each matched pair
    bar1,bar2 = bar.split('-')                                  # split barcode into 1st and 2nd indices
    val=match_pairs[bar]
    percentage=np.round(100*(val/total_reads),2)
    f.write(f"\n{bar1}\t{bar2}\t{val}\t{percentage}")
for bar in ihopped_pairs:                                       # write stats for each i-hopped pair
    bar1,bar2 = bar.split('-')
    val=ihopped_pairs[bar]
    percentage=np.round(100*(val/total_reads),2)
    f.write(f"\n{bar1}\t{bar2}\t{val}\t{percentage}")
f.close()

### Print total counts for each category
print(f"Total Reads: {total_reads}")
print(f"Total Matched Reads: {matches} ({np.round(100*(matches/total_reads),2)}%)")
print(f"Total Index-Hopped Reads: {ihops} ({np.round(100*(ihops/total_reads),2)}%)")
print(f"Total Unknown Reads: {unknowns} ({np.round(100*(unknowns/total_reads),2)}%)")
print(f"Total Reads Below Cutoff: {lowq} ({np.round(100*(lowq/total_reads),2)}%)")
