#!/usr/bin/env python

import bioinfo  as bi
import gzip
import matplotlib.pyplot as plt
import numpy as np

path = "/projects/bgmp/shared/2017_sequencing/"
fnames = ['1294_S1_L008_R1_001.fastq.gz',
'1294_S1_L008_R2_001.fastq.gz',
'1294_S1_L008_R3_001.fastq.gz',
'1294_S1_L008_R4_001.fastq.gz']

# fnames = ['R1_trimmed.fq.gz',
#           'R2_trimmed.fq.gz',
#           'R3_trimmed.fq.gz',
#           'R4_trimmed.fq.gz']

q_dists = {'read1':[[0 for i in range(42)] for j in range(101)],
           'read2':[[0 for i in range(42)] for j in range(101)], 
           'index1':[[0 for i in range(42)] for j in range(8)],
           'index2':[[0 for i in range(42)] for j in range(8)]}
q_sums = {'read1':[0 for i in range(101)],
           'read2':[0 for i  in range(101)], 
           'index1':[0 for i  in range(8)],
           'index2':[0 for i  in range(8)]}
recordctr = 0
with gzip.open(path+fnames[0],'rt') as R1,\
gzip.open(path+fnames[1],'rt') as R2,\
gzip.open(path+fnames[2],'rt') as R3,\
gzip.open(path+fnames[3],'rt') as R4:
    prev = ''
    for l in R1:
        line = l.strip()
        if prev == '+':
            for i in range(101):
                qscore = bi.convert_phred(line[i])
                q_sums['read1'][i]+=qscore
                q_dists['read1'][i][qscore]+=1
            recordctr+=1
        prev = line
    for l in R2:
        line = l.strip()
        if prev == '+':
            for i in range(8):
                qscore = bi.convert_phred(line[i])
                q_sums['index1'][i]+=qscore
                q_dists['index1'][i][qscore]+=1
        prev = line
    for l in R3:
        line = l.strip()
        if prev == '+':
            for i in range(8):
                qscore = bi.convert_phred(line[i])
                q_sums['index2'][i]+=qscore
                q_dists['index2'][i][qscore]+=1
        prev = line
    for l in R4:
        line = l.strip()
        if prev == '+':
            for i in range(101):
                qscore = bi.convert_phred(line[i])
                q_sums['read2'][i]+=qscore
                q_dists['read2'][i][qscore]+=1
        prev = line

q_means = {}
for key in q_sums:
    q_means[key] = []
    for i in range(len(q_sums[key])):
        q_means[key].append(q_sums[key][i]/recordctr)
     
    plt.plot(q_means[key])
    plt.title("Mean Q-Scores by Base Pair #")
    plt.xlabel("Base Pair #")
    plt.ylim(0,50)
    plt.ylabel("Mean Q-Score")
    plt.savefig(f"{key}_means")
    plt.clf()
