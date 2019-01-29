
# coding: utf-8

# In[110]:


#fasta_filename = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short.fasta"


# In[109]:


class Protein(object):
    
    def __init__(self,identifier,sequence):
        self.identifier = str(identifier)
        self.sequence = str(sequence)
        
    def get_identifier(self):
        return self.identifier
    
    def get_sequence(self):
        return self.sequence
    
    def get_mw(self):
        aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}
        return float(sum(map(lambda aa: aminoacid_mw[aa] if aa in aminoacid_mw else 0, self.sequence)))
    
    def has_subsequence(self,subsequence_protein):
        if subsequence_protein.get_sequence() in self.sequence: 
            return True
        else: 
            return False
        
    def get_length(self):
        return int(len(self.get_sequence()))


# In[121]:


def FASTA_iterator(fasta_filename): 
    file = open(fasta_filename, "r")
    sequence =""
    identifier = file.readline().replace(">", "").strip()
    for line in file: 
        if line.startswith(">"):
            yield Protein(identifier, sequence)
            sequence = "" 
            identifier = line.replace(">", "").strip()
        else: 
            sequence+= line.strip()
    file.close()
    yield Protein(identifier, sequence)

