
# coding: utf-8

# In[270]:


def count_sequences_by_residue_threshold(filename,residue,threshold=0.03):
    """
    returns an integer corresponding to the total number of protein sequences having 
    a relative frequency higher or equal than a given threshold for a given residue.     
    """
    count=0
    length_of_current_protein = 0
    frequency_of_residue = 0   
    
    with open(filename, "r") as fd:
        #Skip the first line, not to have an error for the values of the frequency and total sum are still equal 0
        next(fd)
        for line in fd:
            #New Protein found, save data for the old one!
            if(line.startswith(">")):
                    relative_frequency = frequency_of_residue/length_of_current_protein
                    if( relative_frequency >= threshold ):
                        count+=1
                    # Re-initialize             
                    length_of_current_protein = 0
                    frequency_of_residue = 0
                    continue
                    
            length_of_current_protein += len(line.strip())
            frequency_of_residue += line.count(residue)   
        fd.close()
    
    # Catch last case 
    relative_frequency = frequency_of_residue/length_of_current_protein
    if( relative_frequency > threshold ):
         count+=1
    
    return count
    


# In[219]:


# ----- SOLUTION WITH BIOPYTHON -------

#!pip install biopython
# from Bio import SeqIO
# def count_sequences_by_residue_threshold(filename,residue,threshold=0.03):
#     count = 0 
#     with open(filename, mode='r') as handle:
#         for record in SeqIO.parse(handle, 'fasta'):
#             if(record.seq.count(residue)/len(record.seq) > threshold) : count+=1
#     return count


# In[303]:


from collections import OrderedDict
def create_string(identifier ,residues, first_n,last_m):
    first = residues[:first_n]
    last = residues[-last_m:]    
    frequencies=""
    header=""
    for charachter in OrderedDict.fromkeys(first+last):
        frequencies+= header+charachter+":"+str(residues.count(charachter))
        header=","
    output_string=identifier+"\t"+first +"\t"+last+"\t"+frequencies+"\n"
    return output_string

def print_sequence_tails(filename,output_filename,first_n=10,last_m=10):
    """
    Given a protein FASTA file (filename), save on a file (output_filename) 
    the protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency 
    in the protein of each of the first N-aminoacids and the last M-aminoacids. The first three fields 
    must be separated by a tabulator, and the absolute frequency of the residues must have the format 
    RESIDUE:frequency and must be separated by comma. One protein by line. The first line must be a line with
    a summary formatted as follows 
    """
    f = open(output_filename, "w")
    # First line summary 
    n_proteins=1
    identifier=""
    residues=""
    with open(filename, "r") as fd: 
        identifier=fd.readline().replace(">","").strip()
        for line in fd:
            #New Protein found, print the old one!
            if(line.startswith(">")):
                n_proteins+=1
                #Save 
                f.write( create_string(identifier, residues ,first_n,last_m))
                # Re-initialize
                identifier=line.replace(">","").strip()
                residues=""
            else:
                residues+=line.strip()
                
    f.write( create_string(identifier, residues,first_n,last_m))
    f.close()
    
    # Prepend the summary line 
    test_string="#The file "+filename+" contains "+str(n_proteins)+" proteins. Here we show the code of the protein, the first "+str(first_n)+" aminoacids of each protein and the last "+str(last_m)+" aminoacids.\n"
    with open(output_filename, 'r') as original: data = original.read()
    with open(output_filename, 'w') as modified: modified.write(test_string + data)


# In[304]:


#file="/Users/luisasantus/Desktop/PYT/PythonCourse/data/example.fasta"
#count_sequences_by_residue_threshold(file,"C",threshold=0.03)


# In[306]:


#my_file="/Users/luisasantus/Desktop/PYT/PythonCourse/Ipython/example.fasta"
#my_file_out="/Users/luisasantus/Desktop/PYT/PythonCourse/Ipython/example_res.fasta"
#print_sequence_tails(my_file,my_file_out,3,5 )


# In[284]:


def print_sequence_tails_biopython(filename,output_filename,first_n=10,last_m=10):
    """
    Given a protein FASTA file (filename), save on a file (output_filename) 
    the protein identifier, the first N-aminoacids, the last M-aminoacids and the absolute frequency 
    in the protein of each of the first N-aminoacids and the last M-aminoacids. The first three fields 
    must be separated by a tabulator, and the absolute frequency of the residues must have the format 
    RESIDUE:frequency and must be separated by comma. One protein by line. The first line must be a line with
    a summary formatted as follows 
    """
    f = open(output_filename, "w")
    identifier=""
    residues=""
    
    with open(filename, "r") as fd: 
        identifier=fd.readline().replace(">","").strip()
        for record in SeqIO.parse(fd, "fasta"):
                f.write( create_string(record.identifier, record.seq  ,first_n,last_m))
    f.write( create_string(record.identifier, record.seq  ,first_n,last_m))
    f.close()

