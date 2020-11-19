#!/usr/bin/env python
# coding: utf-8

# # Title
# 
# ## Purpose
# State the purpose of the notebook.
# 
# ## Methodology
# Quickly describe assumptions and processing steps.
# 
# ## WIP - improvements
# Use this section only if the notebook is not final.
# 
# Notable TODOs:
# - todo 1;
# - todo 2;
# - todo 3.

# # Setup
# 
# ## Library import
# We import all the required Python libraries

# In[67]:


# Data manipulation
import pandas as pd
import numpy as np

# Options for pandas
pd.options.display.max_columns = 50
pd.options.display.max_rows = 30

# Visualizations
import matplotlib.pyplot as plt

# Include local library paths
import sys
import os 
import re


# 
# # Data import
# We retrieve all the required data for the analysis.

# In[116]:


data=pd.read_csv("./algorithms_computational_biology/T6/expressionExampleGSE121594.tsv",skiprows=0,delimiter="\t")

data = data.set_index("ID")

data = data.loc[~(data==0).all(axis=1)]

coexp_mat = data.T.corr()
df = coexp_mat
np.fill_diagonal(df.values, 0)


# In[117]:


def generate_interaction_list( coexp_mat, min_value = 0.7 ):
    with open("data_coex.txt","w") as data_cyt:
        data_adj_mat = coexp_mat.values
        for i in range(len(data_adj_mat)):
            for j in range(i + i,len(data_adj_mat)):
                if( abs(data_adj_mat[i][j])> min_value and list(coexp_mat.index)[i] != list(coexp_mat.index)[j]):
                    print(list(coexp_mat.index)[i],"\t",list(coexp_mat.index)[j],"\t",data_adj_mat[i][j])
                    inter=list(coexp_mat.index)[i]+"    "+list(coexp_mat.index)[j]+"    "+str(data_adj_mat[i][j])+"\n"
                    data_cyt.write(inter)
                    #data_cyt.write((list(coexp_mat.index)[i],"\t",list(coexp_mat.index)[j],"\t",data_adj_mat[i][j],"\n"))


# In[118]:


#generate_interaction_list(coexp_mat)


# In[105]:


df[df> 0.8 ] = 1
df[df< 0.8 ] = 0


# In[106]:

#print("df",df)


# In[108]:


degrees=df.sum(axis=1).values
#print(degrees)


# In[111]:



from collections import Counter
degres_allprot=dict(Counter(degrees))
"""
plt.figure()
plt.title("Degree Distribution")
plt.xlabel("degree")
plt.ylabel("# of nodes")
plt.bar(degres_allprot.keys(),degres_allprot.values(),)
plt.show()
"""

# In[119]:

"""
edges = np.triu(df.values).sum()#-np.trace(df.values)
print("edges")

total_density = 2*len(df)/(edges*(edges-1))
print("total density")
print(total_density)
print(df.values)
"""

def adjacency_to_dict( adj_matrix ):
    dic = {}
    for i in range(len(adj_matrix)):
        for j in range(len(adj_matrix[i])):
            if adj_matrix[i][j] == 1.0:
                if i in dic:
                    dic[i][j] = 1
                else:
                    dic[i] = {j:1}
    return  dic

def spectrum():

    adj_dict = adjacency_to_dict( df.to_numpy() )
    degree_dict = {}
    print("started!")
    for node_key in adj_dict:
        matches = 0
        for edge_key in adj_dict[ node_key ]: #see which edges share the same edge
            for candidate in adj_dict[ edge_key ]:
                if candidate in adj_dict[ node_key ]:
                    matches += 1

        degree = len( list( adj_dict[ node_key ] ) )
        if degree == 1:
            grouping_coefficient = 0
        else:
            grouping_coefficient = ( matches / 2)  / ( (degree * (degree - 1)) / 2)  # matches / posible matches, same as local_density
        if degree in degree_dict :
            degree_dict[ degree ][ 'list' ].append( grouping_coefficient )
            degree_dict[degree][ "sum" ] += grouping_coefficient
        else:
            degree_dict[ degree ] = {
                "list" : [ grouping_coefficient ],
                "sum" : grouping_coefficient
            }
        print(node_key)

    averages = []
    for degree in degree_dict:
        averages.append( degree_dict[ degree ][ "sum" ] / len( degree_dict[ degree ][ "list" ] ) )

    plt.figure()
    plt.title( "average local density vs degree" )
    plt.xlabel( "Degree" )
    plt.ylabel( "average local density" )
    plt.bar( list( degree_dict.keys() ), averages,  )
    plt.show()
    return degree_dict



    #spectrum( degres_allprot.keys(), [ node[2] for node in file.to_numpy() ] )
spectrum()



# In[113]:


#import networkx as nx


# In[114]:


#G = nx.from_pandas_adjacency(df)






# In[115]:


#print("info",nx.info(G))


# In[93]:


#print("centratily",nx.degree_centrality(G))


# In[94]:


#print("clustering",nx.clustering(G))


# In[1]:




# In[23]:


#a=data.corr().to_numpy()


# In[26]:


#print("a",a)


# # Data processing
# Put here the core of the notebook. Feel free di further split this section into subsections.

# In[ ]:





# In[ ]:




