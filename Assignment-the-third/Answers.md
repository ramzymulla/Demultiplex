Per-barcode counts/percentages: [stats.txt](./dplexer_out/stats.txt)

Summary counts/percentages: [summary.txt](./dplexer_out/summary.txt)
```
Total Reads: 363246735
Total Matched Reads: 304980270 (83.96%)
Total Index-Hopped Reads: 517612 (0.14%)
Total Unknown Reads: 57748853 (15.9%)
Total Reads Below Cutoff: 26964891 (7.42%)
```

Confirmation of correct read count:
```
$ zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz | grep -c "^+"
363246735
```

Q-score Cutoff Rationale:
```
I chose to use a cutoff of Q30 because, from inspetion of the Q-score distributions, the lowest average Q-score (bp #0) was just over Q30, with the aveage over all bp #s being around Q35. This indicated that Q30 would be a solid option for capturing a large proportion (if not all) of the low-quality reads without cutting too far into the pool of usable reads. Further, since Q30 means there is a 1/1000 chance of a specific base being incorrect, an average >=Q30 suggests an extremely low probability of there being any errors in a valid 8-bp index, especially since multiple SNPs are required to bridge the edit-distance between any two valid indices. Additionally, Q30 is widely considered to be the benchmark for read quality in next-gen sequencing (per Illumina website). 
```
