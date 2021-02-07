import numpy as np
import src
from src import Blosum50
from src import GeneticCode
from src import Trie
from src import OptimalAlignmentProtein

class BasicLocalAlignmentSearchTool(object):
    """
    Tool for searching local alignments between two sequences

    Parameters:
    g : str - sequence to be aligned
    t : str - sequence where local alignments g will be found
    th : float - threshold of detection for perfect matchings
    
    """

    def __init__(self, g, t, th):
        self.g = g
        self.t = t
        self.th = th

    def perfect_matches(self, g, t, th, k, open_gap = 10, increase_gap = 1):
        """
        Computes all the indices that correspond to beginning of perfect matches between an element of  Sg  and a subword of t.

        Parameters:
        g : str - sequence to be aligned
        t : str - sequence where local alignments g will be found
        th : float - threshold of detection for perfect matchings
        open_gap : int - penalty of gap opening
        increase_gap : int - penalty of gap increasing

        Returns:

        S_g : Trie, index_list : list - perfect matches between an element of  Sg  and a subword of t, the indices that correspond to beginning of the perfect matches 
        """
        W_g = Trie.Trie()
        S_g = Trie.Trie()
        index_list = []
        for i in range(len(g)-k+1): #cycle over all k-subwords of g
          k_word = g[i:i+k]
          W_g.insert(k_word)
    
        for i in range(len(t)-k+1): #cycle over all k-subwords of t
          k_word_t = t[i:i+k]
          for word in W_g.all_words_iterator(W_g.root):
            self_score = OptimalAlignmentProtein.OptimalAlignmentProtein.optimal_alignment_refined(word, word, open_gap, increase_gap)
            if OptimalAlignmentProtein.OptimalAlignmentProtein.optimal_alignment_refined(k_word_t, word, open_gap, increase_gap) >= th*self_score:
              index_list.append(i)
              S_g.insert(k_word_t)
    
        return S_g, index_list

    def max_alignment(self, g, t, th,th_l = 0.1, k = 4, open_gap = 10, increase_gap = 1):
        """
        Computes all the local alignments with sufficiently high scores.

        Parameters:
        g : str - sequence to be aligned
        t : str - sequence where local alignments g will be found
        th : float - threshold of detection for perfect matchings
        open_gap : int - penalty of gap opening
        increase_gap : int - penalty of gap increasing
        th_l : float - threshold for alignment accepting

        Returns:

        alignment indices : list - list of [i_left,i_right],[j_left,j_right] where i_left - begining of local alignment on g
                                                                                   i_right - end of local alignment on g
                                                                                   j_left - begining of local alignment on t
                                                                                   j_right - end of local alignment on t
        """
        S_g,index_list = self.perfect_matches(g,t,th,k)

        g_score = OptimalAlignmentProtein.OptimalAlignmentProtein.optimal_alignment_refined(g, g, open_gap, increase_gap)
        
        scores = Blosum50.Blosum50()

        alignment_indices = []
        #for each perfect match between Sg and a subword of t
        l=0
        while (l < len(index_list)):
            t_i = index_list[l]
            k_word_t = t[t_i:t_i+k]
            # find an index in g which corresponds to perfect match
            for g_i in range(len(g)-k+1): #cycle over all k-subwords of g
                k_word_g = g[g_i:g_i+k]
                self_score = OptimalAlignmentProtein.OptimalAlignmentProtein.optimal_alignment_refined(k_word_g, k_word_g,open_gap, increase_gap)
                alignment_score = OptimalAlignmentProtein.OptimalAlignmentProtein.optimal_alignment_refined(k_word_g, k_word_t,open_gap, increase_gap)
                if alignment_score >= th*self_score:
                    #when we found an index we try to extend alignment on the left
                    i_left = g_i - 1
                    j_left = t_i - 1        
        
                    while i_left > 0 and j_left > 0:
                        adding_score = scores.get_score(t[j_left], g[i_left])
                        if adding_score > 0:
                            alignment_score += adding_score
                            i_left -= 1
                            j_left -= 1
                        else:
                            break
        
                    #we need to save indices so they must be inside the strings
                    if i_left < 0 or j_left < 0:
                        i_left += 1
                        j_left += 1

                    #and on the right
                    i_right = g_i + k
                    j_right = t_i + k
  
                    while i_right < len(g) and j_right < len(t):
                        adding_score = scores.get_score(t[j_right], g[i_right])
                        if adding_score > 0:
                            alignment_score += adding_score
                            i_right += 1
                            j_right += 1
                        else:
                            break
          
        
                    #if extended alignment was good enough:
                    if alignment_score > th_l*g_score:
                        alignment_indices.append([[i_left,i_right],[j_left,j_right]])
                        g_i = i_right
                        while index_list[l] < j_right and l<len(index_list)-1:
                            l+=1
            l+=1

        return alignment_indices
