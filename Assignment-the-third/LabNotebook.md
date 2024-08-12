
# \<Project Name> Demultiplex: Assingment the Third Lab Notebook
## Date: 2024-08-02

### Objective
- Write up demultiplex.py script
- Run demultiplex.py on real data files
### Methods
***All work done in base environment, which has python 3.12.3, with numpy and matplotlib installed***
- Got started working on demultiplex.py script
- Initially wrote it to iterate over file lines in the main loop, but switched it to use get_record() to extract the entire record each time (accepts the open file handle as param)
- Ported over my working script to use for generating the unit tests in assignment the first
    - provided opportunity to test the basic logic
    - used boolean variables to just take the first match, index-hopped, unknown (one with "N" and one low q-score`)
### Next Steps
- Finish working on demultiplex.py
- Run script on test files
---

## Date: 2024-08-05

### Methods
- Transferred over the function changes from test.py and tweaked it to work for the actual demultiplex algorithm
- Began testing demultiplex.py on unit tests and was confused to find it was failling to correctly calculate the reverse complement, and could not identify the index-hopped read
    - After trying to failing to fix the issue in the python script, I inspected the unit test files and found that I incorrectly identified and SNP as an index-hop in the test.py script, in addition to writing the reverse complemented R3 sequence instead of the original
### Next Steps
- Fix unit tests
- Run demultiplex.py
---

## Date: 2024-08-06

### Methods
- Fixed the unit tests and pushed the changes to github
- Ran demultiplex.py on the new unit tests and got the correct results
    - Only issue was in file naming of matched outputs (easy fix)
- Changed sum printing to give a tab-delimited output rather than plain english
- Made bash script ```demultiplexer.sh``` to for doing my sbatches
- Submitted real run to sbatch, with outputs directed to ```./Assignment-the-third/dplexer_outputs/```
### Next Steps
- Clean up file directory
- Make Answers.md
- Push final repository
---

## Date: 2024-08-08~12
NOTE: My code worked without issue on the first run, but I did some extra optimizing for nerd reasons, so I am consolidating those changes in one block for August 8 through 12
### Methods
- sbatch output fit expected outcomes (total read counts match total records in og files)
    - read counts/percentages roughly match up to those of my peers
```
Command being timed: "python demultiplex.py -q 30 -1 1294_S1_L008_R1_001.fastq.gz -2 1294_S1_L008_R2_001.fastq.gz -3 1294_S1_L008_R3_001.fastq.gz -4 1294_S1_L008_R4_001.fastq.gz -i ../Assignment-the-first/indexes.txt -p /projects/bgmp/shared/2017_sequencing/ -o ./dplexer_outputs/"
	User time (seconds): 8295.45
	System time (seconds): 54.57
	Percent of CPU this job got: 92%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 2:30:21
```
```
Total Matched Reads: 303645222 (83.59%)
Total Index-Hopped Reads: 513040 (0.14%)
Total Unknown Reads: 59088473 (16.27%)
Total Reads Below Cutoff: 28304511 (7.79%)
```
- I got very annoyed with how long my code took to run, so I began trying to optimize my code by:
    - Creating a revcomp dictionary to more quickly convert known indices, only using the actual revcomp() function for unknown indices (ie. not in indexes.txt)
    - Switching the order of if statements to use "in [dict/set]" as much as possible, only calculating qual_scores when absolutely necessary
- My code was still taking very long to run. After 4+ runs with various changes, I realized that my check_qscore()function was extremely bad. Fixing it got my run time down to under 60 minutes. 

Before
```py
def check_qscore(qscores: str,cutoff:int)->bool:
    mean = np.mean([bi.qual_score(i) for i in qscores])
    return bool(mean > cutoff)
```
After
```py
def check_qscore(qscores: str,cutoff:int,enc=33)->bool:
    num = 0                                     # initialize accumulator
    for i in qscores:
        num += bi.convert_phred(i,offset=enc)   # convert each character and add to num
    return num/len(qscores) >= cutoff
```
```
	Command being timed: "python demultiplex.py -q 30 -1 1294_S1_L008_R1_001.fastq.gz -2 1294_S1_L008_R2_001.fastq.gz -3 1294_S1_L008_R3_001.fastq.gz -4 1294_S1_L008_R4_001.fastq.gz -i ../Assignment-the-first/indexes.txt -p /projects/bgmp/shared/2017_sequencing/ -o ./dplexer_out/"
	User time (seconds): 2853.00
	System time (seconds): 54.13
	Percent of CPU this job got: 91%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 53:02.28
```
- Along the way, I also decided to make some additional changes, which included switch check_qscores() to ">=" instead of strictly ">". I think this fits more with the spirit of the term "cutoff". 
    - Also switched so it writes the per-base counts to stats.txt instead of just printing it
- Predictably, this impacted the final output:
```
Before:
Total Matched Reads: 303645222 (83.59%)
Total Index-Hopped Reads: 513040 (0.14%)
Total Unknown Reads: 59088473 (16.27%)
Total Reads Below Cutoff: 28304511 (7.79%)

After:
Total Reads: 363246735
Total Matched Reads: 304980270 (83.96%)
Total Index-Hopped Reads: 517612 (0.14%)
Total Unknown Reads: 57748853 (15.9%)
Total Reads Below Cutoff: 26964891 (7.42%)
```

- Once I was happy with the runtime, I wrote up my Answers.md, added in-line comments to demultiplex.py, cleaned up the directory, and pushed the final repository to github