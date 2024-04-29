#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import pandas as pd
import torch
from torch_geometric.data import Data
from Utils import compute_KNN_graph

def mydatalist (subjectlist,labellist):
    datalist=[]
    for subject, label in zip(subjectlist,labellist):
        #print(subject,label)
        networkFile=('/booboo_workspace/mlewis/UKB_HCP_Neocortex/SC/raw/'+subject+'_SC.txt')
        #print(networkFile)
        x=torch.tensor(np.loadtxt(networkFile,delimiter=','),dtype=torch.float)
        adjFile=('/booboo_workspace/mlewis/UKB_HCP_Neocortex/SC/StreamlineCount_1/'+subject+'_SC_sc1.txt')
        adj=np.loadtxt(adjFile,delimiter=',')
        edge_index=torch.tensor(np.nonzero(adj),dtype=torch.long)
        edge_attr=torch.tensor(np.ones(edge_index.shape[1]),dtype=torch.float)
        y=torch.tensor(label,dtype=torch.long)
        datalist.append(Data(x=x,edge_index=edge_index,edge_attr=edge_attr,y=y))
        #print(len(datalist))
    return datalist

