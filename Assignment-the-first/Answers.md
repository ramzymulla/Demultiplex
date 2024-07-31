# Assignment the First

## Part 1
1. Be sure to upload your Python script. Provide a link to it here:

| File name | label | Read length | Phred encoding |
|---|---|---|---|
| 1294_S1_L008_R1_001.fastq.gz |  |  |  |
| 1294_S1_L008_R2_001.fastq.gz |  |  |  |
| 1294_S1_L008_R3_001.fastq.gz |  |  |  |
| 1294_S1_L008_R4_001.fastq.gz |  |  |  |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    2. **YOUR ANSWER HERE**
    3. **YOUR ANSWER HERE**

```
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | sed -n '2~4p' | grep -E -c "N"
3976613
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz | sed -n '2~4p' | grep -E -c "N"
3328051
```
## Part 2
1. Define the problem
```
Problem is to take 4 FASTQ files (2 biological reads, 2 index reads) and do the following: (1) incorporate index reads to the header lines of all biological reads in the format "index1-index2" and (2) write all biological reads (with the new headers) into new FASTQ files according to their indices.
```
2. Describe output
```
The expected output is 52 FASTQ files—2 for each of the 24 valid index pairs, plus 2 for index-hopped reads, and 2 for unknown reads (indices below the Q-score cutoff or containing "N" reads)
```
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).

4. Pseudocode
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
