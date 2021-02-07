class TrieNode(object):
    """
    Trie node class

    Parameters:
    ch : str - value of node
    """
    def __init__(self,ch = None): 
        self.children = [None]*21
        # isEndOfWord is True if node represents the end of the word 
        self.isEndOfWord = False
        self.char = ch

