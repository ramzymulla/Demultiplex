#!/usr/bin/env python

import bioinfo  as bi
import gzip
import numpy as np

path = "/projects/bgmp/shared/2017_sequencing/"
fnames = ['1294_S1_L008_R1_001.fastq.gz',
'1294_S1_L008_R2_001.fastq.gz',
'1294_S1_L008_R3_001.fastq.gz',
'1294_S1_L008_R4_001.fastq.gz']

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

barcodes = get_indices('indexes.txt')
match,unknown,ihop,lowq = 0,0,0,0
# index_pairs={}
known_pairs=[]
used_pairs=[]

for i in barcodes:
    for j in barcodes:
        known_pairs.append(f"{i}-{j}")
        # index_pairs[f"{i}-{j}"] = 0
num=0
t1 = open('./Test-input_FASTQ/testR1.fastq','w')
t2 = open('./Test-input_FASTQ/testR2.fastq','w')
t3 = open('./Test-input_FASTQ/testR3.fastq','w')
t4 = open('./Test-input_FASTQ/testR4.fastq','w')
t5 = open('./Test-output_FASTQ/testRecord1.fastq','w')
t6 = open('./Test-output_FASTQ/testRecord2.fastq','w')
outfiles=[t1,t2,t3,t4,t5,t6]
first=True
num_recs=1
with gzip.open(path+fnames[0],'rt') as R1,\
gzip.open(path+fnames[1],'rt') as R2,\
gzip.open(path+fnames[2],'rt') as R3,\
gzip.open(path+fnames[3],'rt') as R4:
    
    handles = [R1,R2,R3,R4]
    while match<num_recs or ihop<num_recs or unknown<num_recs or lowq<num_recs:
        hr1,sr1,qr1=get_record(R1)
        if hr1==0: break
        hi1,si1,qi1=get_record(R2)
        hi2,si2_rc,qi2_r=get_record(R3)
        hr2,sr2,qr2=get_record(R4)
        
        si2 = revcomp(si2_rc)
        qi2 = qi2_r[::-1]

        bar = f"{si1}-{si2}"
        hn1 = f"{hr1} {bar}"
        hn2 = f"{hr2} {bar}"

        record1 = f"{hn1}\n{sr1}\n+\n{qr1}"
        record2 = f"{hn2}\n{sr2}\n+\n{qr2}"

        if 'N' in bar:
            if unknown<num_recs:
                if first:
                    first=False
                else:
                    for t in outfiles: t.write("\n")
                print(f"unknown:\n{record1}\n{record2}\n\n")
                t1.write(f"{hr1}\n{sr1}\n+\n{qr1}")
                t2.write(f"{hi1}\n{si1}\n+\n{qi1}")
                t3.write(f"{hi2}\n{si2_rc}\n+\n{qi2_r}")
                t4.write(f"{hr2}\n{sr2}\n+\n{qr2}")
                t5.write(f"{record1}")
                t6.write(f"{record2}")
                unknown +=1
        elif bar not in known_pairs:
            if unknown<num_recs:
                if first:
                    first=False
                else:
                    for t in outfiles: t.write("\n")
                print(f"unknown:\n{record1}\n{record2}\n\n")
                t1.write(f"{hr1}\n{sr1}\n+\n{qr1}")
                t2.write(f"{hi1}\n{si1}\n+\n{qi1}")
                t3.write(f"{hi2}\n{si2_rc}\n+\n{qi2_r}")
                t4.write(f"{hr2}\n{sr2}\n+\n{qr2}")
                t5.write(f"{record1}")
                t6.write(f"{record2}")
                unknown +=1
        elif not check_qscore(qi1+qi2,30):
            if lowq < num_recs:
                if first:
                    first=False
                else:
                    for t in outfiles: t.write("\n")
                print(f"lowq:\n{record1}\n{record2}\n\n")
                t1.write(f"{hr1}\n{sr1}\n+\n{qr1}")
                t2.write(f"{hi1}\n{si1}\n+\n{qi1}")
                t3.write(f"{hi2}\n{si2_rc}\n+\n{qi2_r}")
                t4.write(f"{hr2}\n{sr2}\n+\n{qr2}")
                t5.write(f"{record1}")
                t6.write(f"{record2}")
                lowq+=1
        elif si1!=si2:
            if ihop <num_recs:
                if first:
                    first=False
                else:
                    for t in outfiles: t.write("\n")
                print(f"ihop:\n{record1}\n{record2}\n\n")
                t1.write(f"{hr1}\n{sr1}\n+\n{qr1}")
                t2.write(f"{hi1}\n{si1}\n+\n{qi1}")
                t3.write(f"{hi2}\n{si2_rc}\n+\n{qi2_r}")
                t4.write(f"{hr2}\n{sr2}\n+\n{qr2}")
                t5.write(f"{record1}")
                t6.write(f"{record2}")
                ihop+=1
        
        elif match < num_recs:
            if first:
                first=False
            else:
                for t in outfiles: t.write("\n")
            print(f"match:\n{record1}\n{record2}\n\n")
            t1.write(f"{hr1}\n{sr1}\n+\n{qr1}")
            t2.write(f"{hi1}\n{si1}\n+\n{qi1}")
            t3.write(f"{hi2}\n{si2_rc}\n+\n{qi2_r}")
            t4.write(f"{hr2}\n{sr2}\n+\n{qr2}")
            t5.write(f"{record1}")
            t6.write(f"{record2}")
            match +=1

t1.close()
t2.close()
t3.close()
t4.close()
t5.close()
t6.close()