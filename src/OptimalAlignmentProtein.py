import numpy as np
import src
from src import Blosum50
from src import GeneticCode

class OptimalAlignmentProtein(object):
    '''
  Computes one optimal alignment between two sequences of amino acids 
  obtained from two sequences of nucleotides given in an input file using Blosum50

  Parameters:

  file_path: str - path to input file with 2 strings

  Methods:

  optimal_alignment() - returns one optimal alignment between two sequences of amino acids
  optimal_alignment_affine_penalty(self, open_gap, increase_gap) - returns one optimal alignment between two sequences with affine penalty


    '''
    def __init__(self,file_path):
        self.first_seq, self.second_seq = self.input_processing(file_path)
    
        genetic_code = GeneticCode.GeneticCode()
        self.first_protein = genetic_code.get_protein(self.first_seq)
        self.second_protein = genetic_code.get_protein(self.second_seq)
        self.scores = Blosum50.Blosum50()

        self.n = len(self.first_protein)
        self.m = len(self.second_protein)

    def optimal_alignment(self):
        '''
        Computes one optimal alignment between two sequences given in an input file
        
        Returns:

        first_seq_aligned, second_seq_aligned : str - aligned versions of sequences

        '''
        first_seq_aligned = ""
        second_seq_aligned = ""

        self.C = np.zeros((self.n+1,self.m+1))    

        for i in range(1,self.n+1):
          for j in range(1,self.m+1):
            self.C[i,j] = np.max([self.C[i-1,j] + self.scores.get_score(self.first_protein[i-1], '-'),
                                  self.C[i,j-1] + self.scores.get_score('-', self.second_protein[j-1]), 
                                  self.C[i-1,j-1]+self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1])])

        i = self.n
        j = self.m 

        while (i > 0) and (j > 0):
          if self.C[i,j] == self.C[i-1,j-1]+self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1]):
            first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
            second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
            i -= 1
            j -= 1        
          else:
            if self.C[i,j] == self.C[i,j-1] + self.scores.get_score('-', self.second_protein[j-1]):
              second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
              first_seq_aligned = '-' + first_seq_aligned
              j -=1
            else:
              first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
              second_seq_aligned = '-' + second_seq_aligned
              i -= 1
    
        while i > 0:
          first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
          second_seq_aligned = '-' + second_seq_aligned
          i -= 1
    
        while j > 0:
          second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
          first_seq_aligned = '-' + first_seq_aligned
          j -= 1

        return first_seq_aligned, second_seq_aligned

    def optimal_alignment_affine_penalty(self, open_gap, increase_gap):
        '''
        Computes one optimal alignment between two sequences given in an input file with affine penalty

        Parameters:

        open_gap : int - penalty of gap opening
        increase_gap : int - penalty of gap increasing
        
        Returns:

        first_seq_aligned, second_seq_aligned : str - aligned versions of sequences

        '''
        first_seq_aligned = ""
        second_seq_aligned = ""

        self.C = np.zeros((self.n+1,self.m+1))
        self.P = np.zeros((self.n+1,self.m+1))
        self.Q = np.zeros((self.n+1,self.m+1))
    
        for i in range(1,self.n+1):
          for j in range(1,self.m+1):
            self.P[i,j] = np.max([self.C[i-1,j] - open_gap - increase_gap,
                                  self.P[i-1,j] - increase_gap])
            self.Q[i,j] = np.max([self.C[i,j-1] - open_gap - increase_gap,
                                  self.Q[i,j-1] - increase_gap])
            self.C[i,j] = np.max([self.C[i-1,j-1] + self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1]),
                                  self.P[i,j],
                                  self.Q[i,j]])    
        i = self.n
        j = self.m 

        while i > 0 and j > 0:
          if self.C[i,j] == self.C[i-1,j-1] + self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1]):
            first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
            second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
            i -= 1
            j -= 1
          if self.C[i,j] == self.P[i,j]:
            first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
            second_seq_aligned = '-' + second_seq_aligned
            i -= 1
          if self.C[i,j] == self.Q[i,j]:
            second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
            first_seq_aligned = '-' + first_seq_aligned
            j -=1
        
        while i > 0:
          first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
          second_seq_aligned = '-' + second_seq_aligned
          i -= 1
    
        while j > 0:
          second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
          first_seq_aligned = '-' + first_seq_aligned
          j -= 1

        return first_seq_aligned, second_seq_aligned

    def local_optimal_alignment(self, open_gap, increase_gap):
        '''
        Computes one local optimal alignment between two sequences given in an input file 
        based on the Blosum50 matrix and with affine gap penalty.

        Parameters:

        open_gap : int - penalty of gap opening
        increase_gap : int - penalty of gap increasing
        
        Returns:

        first_seq_aligned, second_seq_aligned : str - aligned versions of sequences

        '''
        first_seq_aligned = ""
        second_seq_aligned = ""

        self.C = np.zeros((self.n+1,self.m+1))

        max_index_i, max_index_j = 0,0
        max_C = 0
        for i in range(1,self.n+1):
          for j in range(1,self.m+1):
            gap_on_first = self.C[:i,j] - open_gap - np.arange(i-1,-1,-1)*increase_gap
            gap_on_second = self.C[i,:j] - open_gap - np.arange(j-1,-1,-1)*increase_gap
            self.C[i,j] = np.max([self.C[i-1,j-1]+self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1]),
                                 np.max(gap_on_first),
                                 np.max(gap_on_second),
                                 0])
            if self.C[i,j] > max_C:
              max_C = self.C[i,j]
              max_index_i, max_index_j = i,j

        i = max_index_i
        j = max_index_j
    
        while self.C[i,j] != 0:
          if self.C[i,j] == self.C[i-1,j-1] + self.scores.get_score(self.first_protein[i-1], self.second_protein[j-1]):
            first_seq_aligned = self.first_protein[i-1] + first_seq_aligned
            second_seq_aligned = self.second_protein[j-1] + second_seq_aligned
            i -= 1
            j -= 1        
          else:

            for ii in range(i):
              if self.C[i,j] == self.C[ii,j] - open_gap - increase_gap*(i-ii-1):
                for iii in range(i,ii,-1):      
                  first_seq_aligned = self.first_protein[iii-1] + first_seq_aligned
                  second_seq_aligned = '-' + second_seq_aligned
                i=ii
                break

            for jj in range(j):
              if self.C[i,j] == self.C[i,jj] - open_gap - increase_gap*(j-jj-1):
                for jjj in range(j,jj,-1):      
                  first_seq_aligned = self.first_protein[jjj-1] + first_seq_aligned                   
                  second_seq_aligned = '-' + second_seq_aligned
                j=jj
                break

        return first_seq_aligned, second_seq_aligned


    def input_processing(self, file_path):
      '''
      Process input for tasks

      Parameters:

      file_path: str - path to input file with 2 strings

      Returns:

      first_seq, second_seq - two lists of characters
      '''
      with open(file_path, 'r') as file:
        data = file.readlines()
  
      data[0] = data[0].replace('\n', '')

      first_seq = [char for char in data[0]]  
      second_seq = [char for char in data[1]]  
  
      return first_seq, second_seq

    @staticmethod
    def optimal_alignment_refined(first_protein, second_protein, open_gap, increase_gap):
            '''
            Computes one optimal alignment between two sequences given as parameters

            Parameters:
            first_protein : str - first protein
            second_protein : str - second protein
            open_gap : int - penalty of gap opening
            increase_gap : int - penalty of gap increasing
        
            Returns:

            first_seq_aligned, second_seq_aligned : str - aligned versions of sequences

            '''
            first_seq_aligned = ""
            second_seq_aligned = ""
        
            n = len(first_protein)
            m = len(second_protein)

            C = np.zeros((n+1, m+1))
            P = np.zeros((n+1, m+1))
            Q = np.zeros((n+1, m+1))
        
            scores = Blosum50.Blosum50()
            for i in range(1, n+1):
              for j in range(1, m+1):
                P[i,j] = np.max([C[i-1,j] - open_gap - increase_gap,
                                      P[i-1,j] - increase_gap])
                Q[i,j] = np.max([C[i,j-1] - open_gap - increase_gap,
                                      Q[i,j-1] - increase_gap])
                C[i,j] = np.max([C[i-1,j-1] + scores.get_score(first_protein[i-1], second_protein[j-1]),
                                      P[i,j],
                                      Q[i,j]])    
            
            score = C[n,m]
            return score