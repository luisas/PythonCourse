
# coding: utf-8

# In[538]:


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


# In[562]:


def add_list_of_set_to_frequency_count(input_list):
    frequency= {}
    for input_set in input_list:
        for identifier in input_set:
            if identifier not in frequency: 
                frequency[identifier] = 1 
            else: 
                frequency[identifier] += 1 
    return frequency


# In[557]:


def compare_fasta_file_identifiers(fasta_filenames_list):
    
    #Throw exception if input list is not containing at least one file.
    if len(fasta_filenames_list )< 1 :
        return "ERROR : input list's length must be 1 or higher!"
    
    #Initialize 
    result_dictionary = {}
    list_of_set_identifiers = []
    names = []
    specific = {}
    
    # Read every fasta file and save the identifiers in sets.
    # list_of_set_identifiers contains the identifiers 
    for fasta_filename in fasta_filenames_list:
        set_identifiers = set()
        iterator = FASTA_iterator(fasta_filename)
        for fasta_tuple in iterator: 
            set_identifiers.add(fasta_tuple[0].upper())
        names.append(fasta_filename)
        list_of_set_identifiers.append(set_identifiers)
        specific.setdefault(fasta_filename,set_identifiers) 
    
    # Collect informations about the frequency
    frequency = add_list_of_set_to_frequency_count(list_of_set_identifiers)

    #Initialize values
    intersection = list_of_set_identifiers[0]
    union = list_of_set_identifiers[0]
    # Iterate all sets to compute unions, intersetions and find specific identifiers
    for index in range(0,len(list_of_set_identifiers)):
        #INTERSECT
        intersection = intersection.intersection(list_of_set_identifiers[index])
        #UNION
        union = union.union(list_of_set_identifiers[index])
        #SPECIFIC
        # Only consider the other files
        specific_set = set(list_of_set_identifiers[index])
        for index_2 in range(0,len(list_of_set_identifiers)):
            if index is not index_2: 
                specific_set.difference_update(list_of_set_identifiers[index_2])
            specific[names[index]] = specific_set
                
    result_dictionary["intersection"]= intersection
    result_dictionary["union"]= union
    result_dictionary["frequency"] = frequency
    result_dictionary["specific"] = specific
    
    return result_dictionary
    


# In[376]:


#fasta_filenames_list = ["/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short.fasta", "/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short2.fasta","/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short3.fasta"]


# In[377]:


#fasta_filenames_list = ["/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_short.fasta"]


# In[563]:


#file1 = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/sample_fasta_files/sample_fasta1.fa"
#file2 = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/sample_fasta_files/sample_fasta2.fa"
#file3 = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/sample_fasta_files/uniprot_sprot_short.fasta"
#file_list = [file1,file2,file3]
#compare_fasta_file_identifiers(file_list)

