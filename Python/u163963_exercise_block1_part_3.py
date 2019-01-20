
# coding: utf-8

# In[10]:


import re
def calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions,
                                 output_filename):
       
    # Extract the subsequences
    subsequences = {}
    with open(subsequences_filename, "r") as f:
        for line in f: 
            subsequences.setdefault(line.strip(), 0)
    f.close()

    max_len_sub = len(max(subsequences , key=lambda t:len(t)))
    # Go through fasta file
    sequence =""
    nr_proteins =1 
    with open(fasta_filename, "r") as f:
        next(f)
        for line in f: 
            if line.startswith(">"): 
                for sub in subsequences:
                    if(len(re.findall('(?='+sub+')',sequence)) >= number_of_repetitions ): subsequences[sub] += 1
                nr_proteins +=1
                sequence =""
            else:
                sequence += line.strip()
        f.close()
    
    #Save last protein 
    for sub in subsequences:
        if(sequence.count(sub) > number_of_repetitions ): subsequences[sub] += 1

    # Sort & Writing
    out = open(output_filename, "w")
    out.write('{:<30s}{:>8s}\n'.format('#Number of proteins:',str(nr_proteins)))
    out.write('{:<30s}{:>8s}\n'.format('#Number of subsequences:',str(len(subsequences))))
    out.write("#subsequence proportions:\n")
    for sub_tuple in sorted(subsequences.items(), key=lambda item: item[1] , reverse=True) :
        out.write("{:s}\t{:>10d}\t{:>.4f}\n".format(sub_tuple[0], sub_tuple[1],sub_tuple[1]/nr_proteins))
    out.close()       


# In[11]:


#subsequences_filename = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/subsequences.txt"
#fasta_filename = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/example.fasta"
#output_filename="/Users/luisasantus/Desktop/PYT/PythonCourse/data/example_res3.fasta"
#number_of_repetitions = 5
#calculate_aminoacid_frequencies(fasta_filename,subsequences_filename,number_of_repetitions,
#                                 output_filename)

