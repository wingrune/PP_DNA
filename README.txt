Programming project 
From ADN to formation of proteins : how to align sequences ?
========

The programming project aims to fins optimal alignments between DNA and protein sequences by using dynamic programming.
Algorithms are stored in /src. One can find some our results in /results and the report in /report.


Launch tool
--------

One can find the main script in /bin/main.py. One can change default parameters at the beginig of the script. 
The input files are in /bin/example. The script will run over these files. If you don't want to spend a lot of time, 
please, delete SubSeq_Res(>320).txt from this folder. 

To run script with default parameters please run (your working directory should be the main project directory)

python -m bin.main

Requirements
--------
>300 Mb GPU, Python >=3.6, used packages: numpy, matplolib, os, time

Features
--------

- LongestSubSeq.py - computes the length of one longest commun subsequence between two sequences
  given in an input file
- OptimalAlignment.py - Computes and displays (nicely) one optimal alignment between two sequences
  given in an input file
- OptimalAlignmentProtein.py - Computes one optimal alignment between two sequences of amino acids 
  obtained from two sequences of nucleotides given in an input file using Blosum50
- BasicLocalAlignmentSearchTool.py - Tool for searching local alignments between two sequences

Parameters 
--------

input_folder - folder with input .txt files
result_folder - folder where results will be stored
save_images = False - if True save images for Tasks 4-8 (memory consuming).

###Tasks 5-8###
open_gap = 10 - penalty for opening gap
increase_gap = 1 #penalty for increasing gap

###Tasks 7-8##
g = 'ADCRGHC'  - first string of toy example
t = 'EDADCRGNRADACRGHC'  - second string of toy example
k = 4 - number of characters for perfect matches search
th = 0.9 - upper border to accept a perfect match
th_l = 0.1 - lower border to accept alignment

th_protein = 0.6 - upper border to accept a perfect match in test examples


