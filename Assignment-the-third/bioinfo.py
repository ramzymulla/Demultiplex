#!/usr/bin/env python

# Author: Ramzy Al-Mulla rza@uoregon.edu

# Check out some Python module resources:
#   - https://docs.python.org/3/tutorial/modules.html
#   - https://python101.pythonlibrary.org/chapter36_creating_modules_and_packages.html
#   - and many more: https://www.google.com/search?q=how+to+write+a+python+module

'''This module is a collection of useful bioinformatics functions
written during the Bioinformatics and Genomics Program coursework.
You should update this docstring to reflect what you would like it to say'''

__version__ = "D.1"         # Read way more about versioning here:
                            # https://en.wikipedia.org/wiki/Software_versioning

import math

### consts ###
DNA_bases = set("ATGC")
RNA_bases = set("AUGC")
SAMCOLS = ['qname','flag','rname',
           'pos','mapq','cigar','rnext',
           'pnext','tlen','seq','qual']
### dicts ###
DNA_COMPS = {a:b for a,b in zip(["A","T","G","C","N"],["T","A","C","G","N"])}
CODONS = {
#       xTx           xCx           xAx           xGx
# Txx
    'TTT': 'Phe', 'TCT': 'Ser', 'TAT': 'Tyr', 'TGT': 'Cys',     # TxT
    'TTC': 'Phe', 'TCC': 'Ser', 'TAC': 'Tyr', 'TGC': 'Cys',     # TxC
    'TTA': 'Leu', 'TCA': 'Ser', 'TAA': '---', 'TGA': '---',     # TxA
    'TTG': 'Leu', 'TCG': 'Ser', 'TAG': '---', 'TGG': 'Trp',     # TxG
# Cxx
    'CTT': 'Leu', 'CCT': 'Pro', 'CAT': 'His', 'CGT': 'Arg',     # CxT
    'CTC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',     # CxC
    'CTA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',     # CxA
    'CTG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',     # CxG
# Axx
    'ATT': 'Ile', 'ACT': 'Thr', 'AAT': 'Asn', 'AGT': 'Ser',     # AxT
    'ATC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',     # AxC
    'ATA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',     # AxA
    'ATG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',     # AxG
# Gxx
    'GTT': 'Val', 'GCT': 'Ala', 'GAT': 'Asp', 'GGT': 'Gly',     # GxT
    'GTC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',     # GxC
    'GTA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',     # GxA
    'GTG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'      # GxG
}
### functions ###
def convert_phred(letter: str, offset: int=33) -> int:
    """
    Converts a single character into a phred score
    letter: string of ASCII character phred score
    offset: integer offset value (33 by default)
    """
    return ord(letter)-offset

def qual_score(phred_score: str,enc=33) -> float:
    """takes a string of phred+33 ASCII symbols and returns the average quality score"""
    num = 0
    for i in phred_score: 
        num += convert_phred(i,offset=enc)
    return num/len(phred_score)

def validate_base_seq(seq: str, RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    if RNAflag:
        return set(seq.upper()) <= RNA_bases
    else:
        return set(seq.upper()) <= DNA_bases

def gc_content(seq: str) -> float:
    '''Returns GC content of a DNA or RNA sequence as a decimal between 0 and 1.'''
    gc_count=0
    assert validate_base_seq(seq), ""
    for chr in seq.upper():
        if chr in "GC":
            gc_count += 1
    return gc_count/len(seq)

def calc_median(lst: list) -> float:
    '''Given a sorted list, returns the median value of the list'''
    n = len(lst)-1
    if n%2==1:
        return (lst[n//2]+lst[n//2 +1])/2
    else:
        return float(lst[n//2])
    
def oneline_fasta(input: str, output: str):
    '''
    Writes a new FASTA file with one-line sequences (ie. exactly two lines per record)

    Args:
        input (str): path to FASTA fle
        output (str): path to where the output file should be written
    '''
    with open(input,'r') as f, open(output, 'w') as oof:
        line = f.readline().strip()
        oof.write(line+"\n")
        line = f.readline().strip()
        while line: 
            if line[0]==">":
                oof.write("\n"+line+"\n")
            else:
                oof.write(line)
            line = f.readline().strip()
    pass



### argparser template ###
# import argparse as arp
# def get_args():
#     '''Parses user inputs from command line'''
#     parser = arp.ArgumentParser(description="blah blah blah")                   
#     parser.add_argument("-k", "--kmer_size", help="defines k-mer size",         
#                         required=True)
#     parser.add_argument("-f","--filepath", help="defines path of data file",
#                         required=True)
#     parser.add_argument("-t", "--title", help="base for output filenames",
#                         required=True)
#     return parser.parse_args()
# args = get_args()


### unit tests ###
if __name__ == "__main__":
    # write tests for functions above, Leslie has already populated some tests for convert_phred
    # These tests are run when you execute this file directly (instead of importing it)
    assert convert_phred("I") == 40, "wrong phred score for 'I'"
    assert convert_phred("C") == 34, "wrong phred score for 'C'"
    assert convert_phred("2") == 17, "wrong phred score for '2'"
    assert convert_phred("@") == 31, "wrong phred score for '@'"
    assert convert_phred("$") == 3, "wrong phred score for '$'"
    print("Your convert_phred function is working! Nice job")


    assert calc_median([1,2,3]) == 2
    assert calc_median([5,6,7,8]) == 6.5
    assert calc_median([1,1,1,1,1,1,1,1,100]) == 1
    assert calc_median([7]) == 7
    assert calc_median([50,100]) == 75
    print("Your calc_median function is working! Nice job")

    assert gc_content("GCGCGC") == 1
    assert gc_content("AATTATA") == 0
    assert gc_content("GCATCGAT") == 0.5
    print("Your gc_content function is working! Nice job")

    assert validate_base_seq("GCGCGC") == True
    assert validate_base_seq("GAGTGC",RNAflag=True) == False
    assert validate_base_seq("AAUUAUA") == False
    assert validate_base_seq("AAUUAUA",RNAflag=True) == True
    assert validate_base_seq("GcAtcGAT") == True
    assert validate_base_seq("AATAGAT"), "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True), "Validate base seq does not work on RNA"
    assert validate_base_seq("R is the best!")==False, "Not a DNA string"
    assert validate_base_seq("aatagat"), "Validate base seq does not work on lowercase DNA"
    assert validate_base_seq("aauagau", True), "Validate base seq does not work on lowercase RNA"
    assert validate_base_seq("TTTTtttttTTT")
    print("Your validate_base_seq function is working! Nice job")

    assert qual_score("A") == 32.0, "wrong average phred score for 'A'"
    assert qual_score("AC") == 33.0, "wrong average phred score for 'AC'"
    assert qual_score("@@##") == 16.5, "wrong average phred score for '@@##'"
    assert qual_score("EEEEAAA!") == 30.0, "wrong average phred score for 'EEEEAAA!'"
    assert qual_score("$") == 3.0, "wrong average phred score for '$'"
    print("Your qual_score function is working! Nice job")