
# coding: utf-8

# In[9]:


#pdb_file_path = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/1tim.pdb"
#pdb_file_path = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/6a6q.pdb"
#pdb_file_path = "/Users/luisasantus/Desktop/PYT/PythonCourse/data/4y5r.pdb"


# In[10]:


from statistics import *
import math


# In[11]:


def calc_chain_medium(list_of_sets_of_coordinates):
    min_distances = []
    for index1 in range(0,len(list_of_sets_of_coordinates)):
        for index2 in range(0,len(list_of_sets_of_coordinates)):
            if index1 != index2:
                min_distance= calc_min_distance(list_of_sets_of_coordinates[index1],list_of_sets_of_coordinates[index2])
                min_distances.append(min_distance)
    return mean(min_distances)


# In[12]:



def eucledian_distance(c1,c2):
    #c1 = (x1,y1,z1)
    #c2 = (x2,y2,z2)
    return(math.sqrt(sum([(a - b) ** 2 for a, b in zip(c1, c2)])))
    


# In[13]:


def calc_min_distance(set_1,set_2):
    min_distance = float("inf")
    for triple1 in set_1: 
        for triple2 in set_2:
            calculated_distance = eucledian_distance(triple1,triple2)
            if calculated_distance < min_distance: 
                min_distance = calculated_distance
    return min_distance


# In[14]:



def calculate_pdb_chain_mean_minimum_distances(pdb_file_path):
    file = open(pdb_file_path, "r")
    #TODO: save everythin
    current_nr_res = 0
    current_chain = ""
    current_set_coordinates = set()
    list_of_sets_of_coordinates= []
    dict_res= {}
    for line in file: 
        if line.startswith("ATOM"):
            splitted_line = [line[:6], line[6:11], line[12:16], line[17:20], line[21], line[22:26], line[30:38], line[38:46], line[46:54]]
            chain = splitted_line[4]
            nr_res = splitted_line[5].strip()
            coordinates = (float(splitted_line[6]), float(splitted_line[7]), float(splitted_line[8]))

            # Only true in first case
            if current_chain is "" : current_chain = chain  
            if current_nr_res is 0 : current_nr_res = nr_res  


            # Check if residue changed
            if nr_res != current_nr_res: 
                list_of_sets_of_coordinates.append(current_set_coordinates)
                current_set_coordinates = set()
                current_nr_res = nr_res

            current_set_coordinates.add((coordinates))

            if current_chain is not chain:
                    #TODO: calc min in atoms and calc med
                    medium = calc_chain_medium(list_of_sets_of_coordinates)
                    dict_res[current_chain] = medium
                    list_of_sets_of_coordinates = []
                    current_chain = chain 

    list_of_sets_of_coordinates.append(current_set_coordinates)   
    medium = calc_chain_medium(list_of_sets_of_coordinates)
    dict_res[current_chain] = medium
    return dict_res


# In[15]:


#calculate_pdb_chain_mean_minimum_distances(pdb_file_path)

