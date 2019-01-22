
# coding: utf-8

# In[160]:


#fasta_filename = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short.fasta"
#fasta_filename = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/sample_fasta_files/sample_fasta1.fa"


# In[46]:


def FASTA_iterator(fasta_filename): 
    file = open(fasta_filename, "r")
    sequence =""
    identifier = file.readline().replace(">", "").strip()
    for line in file: 
        if line.startswith(">"):
            yield(identifier, sequence)
            sequence = ""; 
            identifier = line.replace(">", "").strip()
        else: 
            sequence+= line.strip()
    file.close()
    yield(identifier, sequence)


# In[47]:


def get_max_sequence_length_from_FASTA_file(fasta_filename): 
    list = [len(sequence) for identifier, sequence in FASTA_iterator(fasta_filename)]
    return max(list)


# In[49]:


def get_min_sequence_length_from_FASTA_file ( fasta_filename ): 
    list = [len(sequence) for identifier,sequence in FASTA_iterator(fasta_filename)]
    return min(list)


# In[61]:


def get_longest_sequences_from_FASTA_file( fasta_filename ):
    max_length = get_max_sequence_length_from_FASTA_file ( fasta_filename )
    list = [(identifier.upper(),sequence) for identifier,sequence in FASTA_iterator(fasta_filename) if len(sequence) == max_length]
    list.sort(key= lambda x: x[0])
    return list


# In[66]:


def get_shortest_sequences_from_FASTA_file( fasta_filename ): 
    min_length = get_min_sequence_length_from_FASTA_file ( fasta_filename )
    list = [(identifier.upper(),sequence) for identifier,sequence in FASTA_iterator(fasta_filename) if len(sequence) == min_length]
    list.sort(key= lambda x: x[0])
    return list


# In[119]:


aminoacid_mw = {'A': 89.09, 'C': 121.16, 'E': 147.13, 'D': 133.1, 'G': 75.07, 'F': 165.19, 'I': 131.18, 'H': 155.16, 'K': 146.19, 'M': 149.21, 'L': 131.18, 'N': 132.12, 'Q': 146.15, 'P': 115.13, 'S': 105.09, 'R': 174.2, 'T': 119.12, 'W': 204.23, 'V': 117.15, 'Y': 181.19}


# In[173]:


def get_weight_seq(sequence):
    return sum(map(lambda aa: aminoacid_mw[aa] if aa in aminoacid_mw else 0, sequence ))

def get_molecular_weights( fasta_filename ): 
    iterator = FASTA_iterator(fasta_filename)
    return dict(((identifier.upper(), get_weight_seq(sequence)) for identifier,sequence in iterator))


# In[174]:


def get_sequence_with_max_molecular_weight( fasta_filename ):
    dict_weigths = get_molecular_weights( fasta_filename )
    for identifier,sequence in FASTA_iterator(fasta_filename):
        if get_weight_seq(sequence) == max_weight:
            return (identifier,sequence)

