
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

---

