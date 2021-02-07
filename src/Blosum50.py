class Blosum50(object):
    """
    Blosum50 substitution matrix
    
    Methodes:
    get_score(a,b) - returns alignment score of a and b
    
    """
    def __init__(self):
        self.scores = []
        self.scores.append([5])
        self.scores.append([-2, 7 ])
        self.scores.append([-1, -1, 7])
        self.scores.append([-2, -2, 2, 8 ])
        self.scores.append([-1, -4, -2, -4, 13  ])
        self.scores.append([ -1, 1, 0, 0, -3, 7 ])
        self.scores.append([ -1, 0, 0, 2, -3, 2, 6 ])
        self.scores.append([ 0, -3, 0, -1, -3, -2, -3, 8])
        self.scores.append([-2, 0, 1, -1, -3, 1, 0, -2, 10 ])
        self.scores.append([-1, -4, -3, -4, -2, -3, -4, -4, -4, 5])
        self.scores.append([-2, -3, -4, -4, -2, -2, -3, -4, -3, 2, 5])
        self.scores.append([-1, 3, 0, -1, -3, 2, 1, -2, 0, -3, -3, 6])
        self.scores.append([-1, -2, -2, -4, -2, 0, -2, -3, -1, 2, 3, -2, 7])
        self.scores.append([-3, -3, -4, -5, -2, -4, -3, -4, -1, 0, 1, -4, 0, 8])
        self.scores.append([-1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10])
        self.scores.append([ 1, -1, 1, 0, -1, 0, -1, 0, -1, -3, -3, 0, -2, -3, -1, 5])
        self.scores.append([0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 2, 5 ])
        self.scores.append([-3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1, 1, -4, -4, -3, 15])
        self.scores.append([-2, -1, -2, -3, -3, -1, -2, -3, 2, -1, -1, -2, 0, 4, -3, -2, -2, 2, 8 ])
        self.scores.append([0, -3, -3, -4, -1, -3, -3, -4, -4, 4, 1, -3, 1, -1, -3, -2, 0, -3, -1, 5 ])
        self.scores.append([-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, 1])
	
	    # A R N D C Q E G H I L K M F P S T W Y V -
    
        self.char_dict = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C':4, 'Q': 5, 
                          'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L':10, 'K': 11,
                          'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17,
                          'Y': 18, 'V': 19, '-': 20}
    def get_score(self,a,b):
        """
        Computes alignment score of a and b

        Parameters:
        a : str - first amino acide
        b : str - second amino acide

        Returns:
        score: int - alignment score of a and b
        
        """
        int_a = self.char_dict[a]
        int_b = self.char_dict[b]
    
        if int_a >= int_b:
          return self.scores[int_a][int_b]
        else:
          return self.scores[int_b][int_a]

