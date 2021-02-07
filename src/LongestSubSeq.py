import numpy as np

class LongestSubSeq(object):
    '''
  Computes the length of one longest commun subsequence between two sequences
  given in an input file

  Parameters:

  file_path: str - path to input file with 2 strings

  Methods:

  length() - returns the length of one longest commun subsequence between two sequences
  sub_sequence() - returns one longest commun subsequence between two sequences
  
  '''
    def __init__(self,file_path):
        self.first_seq, self.second_seq = self.input_processing(file_path)
        self.n = len(self.first_seq)
        self.m = len(self.second_seq)
        
        #compute matrix C for resolving subproblems
        self.C = np.zeros((self.n+1,self.m+1))

        for i in range(1,self.n+1):
            for j in range(1,self.m+1):
                if self.first_seq[i-1] != self.second_seq[j-1]:
                    self.C[i,j] = np.max([self.C[i-1,j], self.C[i,j-1]])
                else:
                    self.C[i,j] = np.max([self.C[i-1,j], self.C[i,j-1], self.C[i-1,j-1]+1])
  
    def length(self):
        '''
        Computes the length of one longest commun subsequence between two sequences
 
        Returns

        length: int - the length of one longest commun subsequence between two sequences
    
        '''
    
        length = int(self.C[self.n,self.m])
  
        return length

    def sub_sequence(self):
        '''
        Computes one longest commun subsequence between two sequences

        Returns

        LCS: str - one longest commun subsequence between two sequences
    
        '''
        i = self.n
        j = self.m
        LCS = ""
        while (i > 0) and (j > 0):
            if self.first_seq[i-1] == self.second_seq[j-1]:
                LCS = self.first_seq[i-1] + LCS 
                i -= 1
                j -= 1
            else:
                if self.C[i,j] == self.C[i-1,j]:
                    i -= 1
                else:
                    j -=1
        return LCS

    def input_processing(self,file_path):
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

