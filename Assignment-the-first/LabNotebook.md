
# Demultiplex "Assignment: The First"<br>Lab Notebook
## Date: 2024-07-25

### Objective
- Perform initial data exploration in data fastq files
- Create representative .fastq's for doing unit tests
- Write pseudocode for demultiplexing the reads
### Methods
```
Data is located in /projects/bgmp/shared/2017_sequencing/
```
- Using ```$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_RX_001.fastq.gz | head``` on each file: 
    - R1 and R4 are the sequence reads, R2 and R3 are the index reads
    - All files use phred+33 encoding since the a "#" is present in the phred lines (which is only present in phred+33)
- Using ```$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_RX_001.fastq.gz | head -n400 | while read line; do echo ${#line}; done | sort -n | uniq -c``` on each file:
    - R1 & R4 have read lengths of 101 bp
    - R2 & R3 have read lengths of 8 bp
- Calculated number of "N" reads in R2 & R3:
```
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -E -c "N"
3976613
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -E -c "N"
3328051
```
### Next Steps
- Write and run Q-score distribution script
---

## Date: 2024-07-30

### Methods
- Completed and ran Q-score distribution script:
```
Command being timed: "python Qdists.py"
	User time (seconds): 28964.86
	System time (seconds): 17.25
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 8:07:41
```
- Figures were output to ```.../Assignment-the-first```
### Next Steps

---

## Date: YYYY-MM-DD

### Objective

### Methods

### Next Steps

---

