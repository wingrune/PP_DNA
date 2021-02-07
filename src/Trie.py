from src import TrieNode

class Trie(object):
    """
    Trie data structure class for words where letters are amino acids
    
    Methods:

    insert(key) - If not present, inserts key into trie. If the key is prefix of trie node, just marks leaf node 
    search(key) - Search key in the trie 
    all_words_iterator(root,word = [],level=0,alpha_size=21) - static method to create iterator over all words in the trie
    all_words(root,word = [],level=0,alpha_size=21) - static method to loop over trie

    """
    
    def __init__(self): 
        self.root = self.getNode() 
        self.char_dict = {'A': 0, 'R': 1, 'N': 2, 'D': 3, 'C':4, 'Q': 5, 
                      'E': 6, 'G': 7, 'H': 8, 'I': 9, 'L':10, 'K': 11,
                      'M': 12, 'F': 13, 'P': 14, 'S': 15, 'T': 16, 'W': 17,
                      'Y': 18, 'V': 19, '-': 20}
        self.index_dict = {0: 'A', 1: 'R', 2: 'N', 3: 'D', 4: 'C', 5: 'Q', 
                      6: 'E', 7: 'G', 8: 'H', 9: 'I', 10: 'L', 11: 'K',
                      12: 'M', 13: 'F', 14: 'P', 15: 'S', 16: 'T', 17: 'W',
                      18: 'Y', 19: 'V', 20: '-'} 

    def getNode(self,ch=None): 
        """Returns new trie node (initialized to NULLs) """
        return TrieNode.TrieNode(ch) 
  
    def _charToIndex(self,ch): 
         """ private helper function,  converts key current character into index """         
         return self.char_dict[ch]
  
  
    def insert(self,key): 
        """
        If not present, inserts key into trie. If the key is prefix of trie node,  
        just marks leaf node 

        Parameters:

        key : str - key to insert into trie

        """
         
        current_node = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            # if current character is not present 
            if not current_node.children[index]: 
                current_node.children[index] = self.getNode(key[level]) 
            current_node = current_node.children[index] 
  
        # mark last node as leaf 
        current_node.isEndOfWord = True
  
    def search(self, key): 
        """
        Search key in the trie 

        Parameters:

        key : str - key to insert into trie

        Returns:

        Returns true if key presents in trie, else false
        """  
     
        current_node = self.root 
        length = len(key) 
        for level in range(length): 
            index = self._charToIndex(key[level]) 
            if not current_node.children[index]: 
                return False
            current_node = current_node.children[index] 
  
        return current_node != None and current_node.isEndOfWord 
    
    @staticmethod
    def all_words_iterator(root,word = [],level=0,alpha_size=21):
      """
      Static method to create iterator over all words in the trie

      Parameters:
      root : TrieNode - rood of the Trie
      word : str - part of word already reconstructed
      level : int - level of depth in the Trie
      alpha_size : int - number of symbols in the alphabet (21 for amino acids + STOP)

      Returns:

      iteraror over all words in the trie


      """
     
        # If node is leaf node, it indicates end of string

      if root.isEndOfWord:
        yield ''.join(word)
      
      for i in range(alpha_size):
      # if NON NULL child is found 
      # add parent key to str and 
      # call the display function recursively 
      # for child node 
        if (root.children[i]):
          if level < len(word):
            word[level] = root.children[i].char 
          else:
            word.append(root.children[i].char)
          yield from Trie.all_words_iterator(root.children[i],word,level+1)
    
    @staticmethod
    def all_words(root,word = [],level=0,alpha_size=21):
      """
      Static method to create iterator over all words in the trie

      Parameters:
      root : TrieNode - rood of the Trie
      word : str - part of word already reconstructed
      level : int - level of depth in the Trie
      alpha_size : int - number of symbols in the alphabet (21 for amino acids + STOP)

      Returns:

      words (str) stored in the Trie

      """
      # If node is leaf node, it indicates end of string
      if root.isEndOfWord:
        print(''.join(word))
      
      for i in range(alpha_size):
      # if NON NULL child is found 
      # add parent key to str and 
      # call the display function recursively 
      # for child node 
        if (root.children[i]):
          if level < len(word):
            word[level] = root.children[i].char 
          else:
            word.append(root.children[i].char)
          Trie.all_words(root.children[i],word,level+1)

