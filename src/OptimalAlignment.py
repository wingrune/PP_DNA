import numpy as np
import matplotlib.pyplot as plt
class OptimalAlignment(object):
    '''
  Computes and displays (nicely) one optimal alignment between two sequences given in an input file

  Parameters:

  file_path: str - path to input file with 2 strings

  Methods:

  optimal_alignment() - returns one optimal alignment between two sequences
  display_alignment(first_seq_aligned, second_seq_aligned,name_to_save) - displays two aligned sequences
  
  '''
    def __init__(self,file_path):
        self.first_seq, self.second_seq = self.input_processing(file_path)
        self.n = len(self.first_seq)
        self.m = len(self.second_seq)
    

    def optimal_alignment(self):
        '''
        Computes one optimal alignment between two sequences given in an input file
        
        Returns:

        first_seq_aligned, second_seq_aligned : str - aligned versions of sequences

        '''

        first_seq_aligned = ""
        second_seq_aligned = ""

        i = self.n
        j = self.m

        #compute matrix C for resolving subproblems
        self.C = np.zeros((self.n+1,self.m+1))
        for i in range(1,self.n+1):
          for j in range(1,self.m+1):
            if self.first_seq[i-1] != self.second_seq[j-1]:
              self.C[i,j] = np.max([self.C[i-1,j], self.C[i,j-1]])
            else:
              self.C[i,j] = np.max([self.C[i-1,j], self.C[i,j-1], self.C[i-1,j-1]+1])
        # traceback C to reconstruct optimal alignment
        while (i > 0) and (j > 0):
          if self.first_seq[i-1] == self.second_seq[j-1]:
            first_seq_aligned = self.first_seq[i-1] + first_seq_aligned
            second_seq_aligned = self.second_seq[j-1] + second_seq_aligned
            i -= 1
            j -= 1
          else:
            if self.C[i,j] == self.C[i-1,j]:
              first_seq_aligned = self.first_seq[i-1] + first_seq_aligned
              second_seq_aligned = '-' + second_seq_aligned
              i -= 1
            else:
              second_seq_aligned = self.second_seq[j-1] + second_seq_aligned
              first_seq_aligned = '-' + first_seq_aligned
              j -=1
        # strings must be the same length
        while i > 0:
          first_seq_aligned = self.first_seq[i-1] + first_seq_aligned
          second_seq_aligned = '-' + second_seq_aligned
          i -= 1
    
        while j > 0:
          second_seq_aligned = self.second_seq[j-1] + second_seq_aligned
          first_seq_aligned = '-' + first_seq_aligned
          j -= 1

        return first_seq_aligned, second_seq_aligned
    
    @staticmethod
    def display_alignment(first_seq_aligned, second_seq_aligned,name_to_save):
        '''
        Displays two aligned sequences

        Parameters:
        first_seq_aligned : str - first aligned sequence
        second_seq_aligned : str - second aligned sequence
        name_to_save : str - full name of the file where the figure will be stored


        '''
        string_1 = first_seq_aligned
        string_2 = second_seq_aligned
        # Create figure and axes
        fig,ax = plt.subplots(figsize=(len(string_1)//3, 5))

        box_matched = {'facecolor':'blue',
                      'edgecolor': 'blue'} #color of box
        left, right = plt.xlim()  
        dx = right/len(string_1)
        for i in range(len(string_1)):
            if string_1[i] == string_2[i]:
                ax.text(dx*i,0.6,string_1[i], 
                bbox = box_matched ,
                color = 'white',    
                fontsize = 14)
                ax.text(dx*i,0.4,string_2[i],
                bbox = box_matched ,
                color = 'white',    
                fontsize = 14)
            else:
                ax.text(dx*i,0.6,string_1[i],fontsize = 14)
                ax.text(dx*i,0.4,string_2[i],fontsize = 14)

    
        #remove axes
        plt.axis('off')
        plt.savefig(name_to_save)
        

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

