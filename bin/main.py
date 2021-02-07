# import 
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import src
from src import LongestSubSeq
from src import OptimalAlignment
from src import OptimalAlignmentProtein
from src import BasicLocalAlignmentSearchTool
from src import Trie
from src import GeneticCode
##############################################################################
# Parameters

input_folder = './bin/example/' #folder with input .txt files
result_folder = './results/' #folder where results will be stored
save_images = False # if True save images for Tasks 4-8 (memory consuming).

# Tasks 5-8
open_gap = 10 #penalty for opening gap
increase_gap = 1 #penalty for increasing gap

# Tasks 7-8
g = 'ADCRGHC' #first string of toy example
t = 'EDADCRGNRADACRGHC' #second string of toy example
k = 4 # number of characters for perfect matches search
th = 0.9 # upper border to accept a perfect match
th_l = 0.1 # lower border to accept alignment

th_protein = 0.6 # upper border to accept a perfect match in test examples


##############################################################################
# Task 2
print("Task 2: Longest subsequense")

#lists for plot data
list_time = np.array([])
list_nm = np.array([])

#loop over file in input directory
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        start_time = time.time() 
        longest_subseq = LongestSubSeq.LongestSubSeq(input_folder + filename)
        LCS = longest_subseq.sub_sequence()
        stop_time = time.time()
        print(f"The length of one longest subsequense {longest_subseq.length()}")
        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, longest_subseq.n * longest_subseq.m)
  
#create figure
fig = plt.figure(figsize = (16,6))
#plot figure
plt.plot(list_nm, list_time, marker = "*", label="Experiment")
plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")

plt.xlabel(r'$nm$', fontsize = 15)
plt.ylabel('Time of execution', fontsize = 15)

plt.legend(fontsize = 15)

plt.grid()

plt.savefig(result_folder + 'time_vs_nm_LCS.png')

# Task 3

print("\nTask 3 : Graphical representation")
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        longest_subseq = OptimalAlignment.OptimalAlignment(input_folder + filename)
        first_seq_aligned, second_seq_aligned = longest_subseq.optimal_alignment()
        #we don't want to display really long sequences
        if longest_subseq.C[longest_subseq.n,longest_subseq.m] < 500:
            name_to_save = result_folder + filename.split('.')[0] + '.png'
            OptimalAlignment.OptimalAlignment.display_alignment(first_seq_aligned, second_seq_aligned,name_to_save)
            print(name_to_save + ' is saved in the result folder: ' + result_folder)

# Task 4

print("\nTask 4 : Substitution matrices")
#lists for plot data
list_time = np.array([])
list_nm = np.array([])
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        start_time = time.time() 
        opt_alignment = OptimalAlignmentProtein.OptimalAlignmentProtein(input_folder + filename)
        first_seq_aligned, second_seq_aligned = opt_alignment.optimal_alignment()
        stop_time = time.time()
        print('1 sequence aligned:' + first_seq_aligned)
        print('2 sequence aligned:' + second_seq_aligned)

        #we save nice figure of alignment
        if save_images:
            name_to_save = result_folder + filename.split('.')[0] + 'Protein.png'
            OptimalAlignment.OptimalAlignment.display_alignment(first_seq_aligned, second_seq_aligned,name_to_save)
            print(name_to_save + ' is saved in the result folder: ' + result_folder)

        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, opt_alignment.n * opt_alignment.m)

if save_images:  
    #create figure
    fig = plt.figure(figsize = (16,6))
    #plot figure
    plt.plot(list_nm, list_time, marker = "*", label="Experiment")
    plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")

    plt.xlabel(r'$nm$', fontsize = 15)
    plt.ylabel('Time of execution', fontsize = 15)

    plt.legend(fontsize = 15)

    plt.grid()

    plt.savefig(result_folder + 'time_vs_nm_OptimalAlignmentProtein.png')

# Task 5

print("\nTask 5 : Affine penalty")
#lists for plot data
list_time = np.array([])
list_nm = np.array([])
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        start_time = time.time() 
        opt_alignment = OptimalAlignmentProtein.OptimalAlignmentProtein(input_folder + filename)
        first_seq_aligned, second_seq_aligned = opt_alignment.optimal_alignment_affine_penalty(open_gap, increase_gap)
        stop_time = time.time()
        print('1 sequence aligned:' + first_seq_aligned)
        print('2 sequence aligned:' + second_seq_aligned)

        #we save nice figure of alignment
        if save_images:
            name_to_save = result_folder + filename.split('.')[0] + 'ProteinAffinePenalty.png'
            OptimalAlignment.OptimalAlignment.display_alignment(first_seq_aligned, second_seq_aligned,name_to_save)
            print(name_to_save + ' is saved in the result folder: ' + result_folder)

        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, opt_alignment.n * opt_alignment.m)
if save_images:  
    #create figure
    fig = plt.figure(figsize = (16,6))
    #plot figure
    plt.plot(list_nm, list_time, marker = "*", label="Experiment")
    plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")

    plt.xlabel(r'$nm$', fontsize = 15)
    plt.ylabel('Time of execution', fontsize = 15)

    plt.legend(fontsize = 15)

    plt.grid()

    plt.savefig(result_folder + 'time_vs_nm_OptimalAlignmentProteinAffinePenalty.png')

# Task 6

print("\nTask 6 : Local Alignment")
#lists for plot data
list_time = np.array([])
list_nm = np.array([])
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        start_time = time.time() 
        opt_alignment = OptimalAlignmentProtein.OptimalAlignmentProtein(input_folder + filename)
        first_seq_aligned, second_seq_aligned = opt_alignment.local_optimal_alignment(open_gap, increase_gap)
        stop_time = time.time()
        print('1 sequence aligned:' + first_seq_aligned)
        print('2 sequence aligned:' + second_seq_aligned)

        #we save nice figure of alignment
        if save_images:
            name_to_save = result_folder + filename.split('.')[0] + 'ProteinLocalAlignment.png'
            OptimalAlignment.OptimalAlignment.display_alignment(first_seq_aligned, second_seq_aligned,name_to_save)
            print(name_to_save + ' is saved in the result folder: ' + result_folder)

        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, opt_alignment.n * opt_alignment.m *  opt_alignment.m)
  
if save_images:
    #create figure
    fig = plt.figure(figsize = (16,6))
    #plot figure
    plt.plot(list_nm, list_time, marker = "*", label="Experiment")
    plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")

    plt.xlabel(r'$nm^2$', fontsize = 15)
    plt.ylabel('Time of execution', fontsize = 15)

    plt.legend(fontsize = 15)

    plt.grid()

    plt.savefig(result_folder + 'time_vs_nm_OptimalLocalAlignmentProtein.png')

# Task 7

print("\nTask 7 : Perfect Matches")

print("\nToy example:")
opt_alignment = BasicLocalAlignmentSearchTool.BasicLocalAlignmentSearchTool(g,t,th)
S_g, index_list = opt_alignment.perfect_matches(g,t,th,k)
print("\ng: " + g)
print("\nt: " + t)
print("S_g: ")
Trie.Trie.all_words(S_g.root)
print("Indices of perfect matches in t: " + str(index_list))


#lists for plot data
list_time = np.array([])
list_nm = np.array([])
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        #input processing
        with open(input_folder + filename, 'r') as file:
            data = file.readlines()
  
        data[0] = data[0].replace('\n', '')

        first_seq = [char for char in data[0]]  
        second_seq = [char for char in data[1]]  

        genetic_code = GeneticCode.GeneticCode()
        first_protein = ''.join(genetic_code.get_protein(first_seq))
        second_protein = ''.join(genetic_code.get_protein(second_seq))
       
        start_time = time.time() 
        opt_alignment = BasicLocalAlignmentSearchTool.BasicLocalAlignmentSearchTool(first_protein,second_protein,th_protein)
        S_g, index_list = opt_alignment.perfect_matches(first_protein,second_protein,th_protein,k)
       
        stop_time = time.time()
        print("S_g: ")
        Trie.Trie.all_words(S_g.root)
        print("Indices of perfect matches in t: " + str(index_list))


        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, len(first_protein) * len(second_protein))
if save_images:
    #create figure
    fig = plt.figure(figsize = (16,6))
    #plot figure
    plt.plot(list_nm, list_time, marker = "*", label="Experiment")
    plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")

    plt.xlabel(r'$nm$', fontsize = 15)
    plt.ylabel('Time of execution', fontsize = 15)

    plt.legend(fontsize = 15)

    plt.grid()

    plt.savefig(result_folder + 'time_vs_nm_BLAST_Perfect_Matches.png')

# Task 8

print("\nTask 8 : Local Alignment")

print("\nToy example:")
opt_alignment = BasicLocalAlignmentSearchTool.BasicLocalAlignmentSearchTool(g,t,th)
alignment_indices = opt_alignment.max_alignment(g,t,th,th_l, k, open_gap, increase_gap)
print("\ng: " + g)
print("\nt: " + t)

print("Local alignments with sufficiently high scores are :")
for alignment in alignment_indices:
  print(g[alignment[0][0]:alignment[0][1]])
  print(t[alignment[1][0]:alignment[1][1]])


#lists for plot data
list_time = np.array([])
list_nm = np.array([])
for filename in os.listdir(input_folder):
    if filename.endswith(".txt"):
        print(filename + ' is opened.')
        #input processing
        with open(input_folder + filename, 'r') as file:
            data = file.readlines()
  
        data[0] = data[0].replace('\n', '')

        first_seq = [char for char in data[0]]  
        second_seq = [char for char in data[1]]  

        genetic_code = GeneticCode.GeneticCode()
        first_protein = ''.join(genetic_code.get_protein(first_seq))
        second_protein = ''.join(genetic_code.get_protein(second_seq))
       
        start_time = time.time() 
        opt_alignment = BasicLocalAlignmentSearchTool.BasicLocalAlignmentSearchTool(first_protein,second_protein,th_protein)
        alignment_indices = opt_alignment.max_alignment(first_protein,second_protein,th_protein,th_l, k, open_gap, increase_gap)
       
        stop_time = time.time()
        print("Local alignments with sufficiently high scores are :")
        for alignment in alignment_indices:
            print(g[alignment[0][0]:alignment[0][1]])
            print(t[alignment[1][0]:alignment[1][1]])


        list_time = np.append(list_time, stop_time - start_time)
        list_nm = np.append(list_nm, len(first_protein) * len(second_protein) *  len(second_protein))

if save_images:  
    #create figure
    fig = plt.figure(figsize = (16,6))
    #plot figure
    plt.plot(list_nm, list_time, marker = "*", label="Experiment")
    plt.plot(list_nm, (list_time[0]/list_nm[0])*list_nm, label="Theory")


    plt.xlabel(r'$nm^2$', fontsize = 15)
    plt.ylabel('Time of execution', fontsize = 15)

    plt.legend(fontsize = 15)

    plt.grid()

    plt.savefig(result_folder + 'time_vs_nm_BLAST.png')
