class GeneticCode(object):
    """
    Decodes a protein from a DNA string

    Methods:
    get_protein(DNA_string) - returns a decoded protein from a DNA string

    """
    def __init__(self):
        self.genetic_code = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
            'TAC':'Y', 'TAT':'Y', 'TAA':'-', 'TAG':'-', 
            'TGC':'C', 'TGT':'C', 'TGA':'-', 'TGG':'W', 
        }

    def get_protein(self, DNA_string):
        """
        Decodes a protein from a DNA string

        Parameters:
        DNA_string : str - string to be decoded

        Returns:

        protein : list - list of amino acids that form the decoded protein

        """
        protein = [] 
        for i in range(0, len(DNA_string)-2, 3): 
            codon = ""
            for j in range(3):
              codon += DNA_string[i+j] 
        
            protein+= self.genetic_code[codon] 
        return protein 

