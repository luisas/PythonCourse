#!/usr/bin/env python
# coding: utf-8

# In[4]:


### NOTE: could be substituted by an import!

#ALPHABETS
protein_letters = 'ACDEFGHIKLMNPQRSTVWY'
rna_letters = 'GAUC'
dna_letters = 'GATC'


#WEIGHTS

protein_weights = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}

rna_weights = {'A': 363.0, 'C': 339.0, 'U': 340.0, 'G': 379.0}

dna_weights = {'A': 347.0, 'C': 323.0, 'T': 322.0, 'G': 363.0}

#COMPLEMENT

dna_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
rna_complement = {'A': 'U', 'C': 'G', 'G': 'C', 'U': 'A'}


#CODON TABLES

rna_table = {'GUC': 'V', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GUU': 'V', 'AAC': 'N', 'AGG': 'R', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'AGA': 'R', 'AAU': 'N', 'ACU': 'T', 'CAC': 'H', 'GUG': 'V', 'CCG': 'P', 'CCA': 'P', 'AGU': 'S', 'CCC': 'P', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'CCU': 'P', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'GAG': 'E', 'UCC': 'S', 'UAC': 'Y', 'CGU': 'R', 'GAA': 'E', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'UCA': 'S', 'AUG': 'M', 'CUG': 'L', 'AUU': 'I', 'CAU': 'H', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'GAC': 'D', 'GUA': 'V', 'UGC': 'C', 'GCU': 'A', 'UGU': 'C', 'CUC': 'L', 'UUG': 'L', 'UUA': 'L', 'GAU': 'D', 'UUC': 'F'}

rna_table_back = {'A': 'GCU', 'C': 'UGU', None: 'UAA', 'E': 'GAG', 'D': 'GAU', 'G': 'GGU', 'F': 'UUU', 'I': 'AUU', 'H': 'CAU', 'K': 'AAG', 'M': 'AUG', 'L': 'UUG', 'N': 'AAU', 'Q': 'CAG', 'P': 'CCU', 'S': 'UCU', 'R': 'CGU', 'T': 'ACU', 'W': 'UGG', 'V': 'GUU', 'Y': 'UAU'}

rna_stop_codons = ['UAA', 'UAG', 'UGA']
rna_start_codons = ['UUG', 'CUG', 'AUG']



dna_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T', 'ATC': 'I', 'ATA': 'I', 'AGG': 'R', 'CCT': 'P', 'AGC': 'S', 'AGA': 'R', 'ATT': 'I', 'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'GGG': 'G', 'GGA': 'G', 'GGC': 'G', 'TAC': 'Y', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GAG': 'E', 'GTT': 'V', 'GAC': 'D', 'GAA': 'E', 'AAG': 'K', 'AAA': 'K', 'AAC': 'N', 'CTC': 'L', 'CAT': 'H', 'AAT': 'N', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'TGT': 'C', 'TCT': 'S', 'GAT': 'D', 'TTT': 'F', 'TGC': 'C', 'TGG': 'W', 'TTC': 'F', 'TCG': 'S', 'TTA': 'L', 'TTG': 'L', 'TCC': 'S', 'ACC': 'T', 'TCA': 'S', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A'}

dna_table_back = {'A': 'GCT', 'C': 'TGT', None: 'TAA', 'E': 'GAG', 'D': 'GAT', 'G': 'GGT', 'F': 'TTT', 'I': 'ATT', 'H': 'CAT', 'K': 'AAG', 'M': 'ATG', 'L': 'TTG', 'N': 'AAT', 'Q': 'CAG', 'P': 'CCT', 'S': 'TCT', 'R': 'CGT', 'T': 'ACT', 'W': 'TGG', 'V': 'GTT', 'Y': 'TAT'}

dna_stop_codons = ['TAA', 'TAG', 'TGA']
dna_start_codons = ['TTG', 'CTG', 'ATG']


# In[128]:


class Sequence(object):
    alphabet = ''
    mw = {}
    
    def __init__(self, identifier, sequence):
        self.__identifier = str(identifier)
        self.__sequence = sequence
        
        for monomer in self.__sequence: 
            if(monomer not in self.alphabet):
                raise ValueError("Impossible to create instance: %s not possible" %monomer)

    def get_identifier(self):
        return self.__identifier

    def get_sequence(self):
        return self.__sequence

    def get_mw(self):
        return sum(self.mw.setdefault(aa, 0) for aa in self.__sequence)

    def has_subsequence(self, subsequence):
        return subsequence.get_sequence() in self.get_sequence()
    
    def __len__(self): 
        return len(self.__sequence)
    
    def __eq__(self, other):
        return self.get_sequence() == other.get_sequence()
    
    def __ne__(self, other):
        return self.get_sequence() != other.get_sequence()

    def __add__(self, other):
        if self.__class__ != other.__class__:
            raise ValueError("Operation not possible: the Sequences you are trying to add are from different classes: %s and %s" %(str(self.__class__), str(other.__class__)))
        sequence = str(self.get_sequence())+ str(other.get_sequence())
        identifier = str(self.get_identifier())+"+"+str(other.get_identifier())
        return self.__class__(identifier,sequence)
    
    def __getitem__(self, key):
        return self.get_sequence()[key]
    
    def __contains__(self, item):
        return item in self.get_sequence()
    
    def __lt__(self, other):
        return self.get_mw() < other.get_mw()
    
    def __le__(self, other):
        return self.get_mw() <= other.get_mw()
    
    def __gt__(self, other):
        return self.get_mw() > other.get_mw()
    
    def __ge__(self,other):
        return self.get_mw() >= other.get_mw()
    
    def __hash__(self):
        string_to_hash = self.get_identifier()+self.get_sequence()
        return string_to_hash.__hash__()


# In[129]:


class ProteinSequence(Sequence):
    alphabet = protein_letters
    mw = protein_weights


# In[130]:


class NucleotideSequence(Sequence):
    complement = {}
    
    start_codons= []
    stop_codons = []
    codons = {}
    
    def translate(self):
        #Initialize vars
        translation_started = False 
        translated = ""
        index_beginning = 0 
        window = 3
        index_end = index_beginning + window
        s = self.get_sequence()
        
        while index_end <= len(s):
            index_end = index_beginning+window
            current_codon = s[index_beginning:index_end]
            
            if current_codon in self.start_codons :
                #print("Beginning Translation! with :",current_codon)
                translation_started = True

            if translation_started:
                index_beginning+= window
                if current_codon in self.codons: 
                    translated += self.codons[current_codon]      
                if current_codon in self.stop_codons:
                    #print("Translation ended with :",current_codon)
                    return ProteinSequence(self.get_sequence(), translated)
            else: 
                index_beginning+= 1
        return ProteinSequence(self.get_sequence(), translated)
        


# In[131]:


class DNASequence(NucleotideSequence):
    alphabet = dna_letters
    mw = dna_weights
    
    start_codons= dna_start_codons
    stop_codons = dna_stop_codons
    codons = dna_table
    
    def transcribe(self): 
        return RNASequence(self.get_identifier(),self.get_sequence().replace("T","U"))


# In[132]:


class RNASequence(NucleotideSequence):
    alphabet = rna_letters
    mw = rna_weights
    
    start_codons= rna_start_codons
    stop_codons = rna_stop_codons
    codons = rna_table
    
    def reverse_transcribe(self):
        return DNASequence(self.get_identifier(),self.get_sequence().replace("U","T"))


# In[149]:


# a = ProteinSequence("identifier1", "AAA")
# t = ProteinSequence("identifier1", "AAA")
# b = DNASequence("identifier2", "TAGATTGTAGTTTAGGAAAATAGAAAAAAAAAAAAAATAG")
# c = DNASequence("identifier3", "AAAATTTTTTTTTCC")
# d = RNASequence("identifier4", "UUAA")
# f = DNASequence("identifier5", "AAAATTCC")

