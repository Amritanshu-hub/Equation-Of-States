#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import csv


# In[2]:


T = 373  # Kelvin
v = 0.0005 # m3
mole = 0.5 
R = 8.314 # J/mol.k
Vm = v/mole


# In[3]:


yi = np.array([0.05, 0.1, 0.05,0.05, 0.1,  0.1, 0.1, 0.05, 0.15, 0.1, 0.1, 0.05])
Tci = np.array([190.6,305.4,370,408.1,425.1,460.8,469.6,507.4,540.6,568.8,595,645.5])
Pci = np.array([4.6,4.88,4.3,3.65,3.8,3.29,3.37,3.3,2.74,2.49,2.3,2.39])
w = np.array([0.0115,0.099,0.153,0.183,0.199,0.227,0.251,0.299,0.349,0.398,0.443,0.484])


# In[4]:


bi = 0.08664*R*Tci/(Pci*pow(10,6))


# In[5]:


Tr = 373/Tci


# In[7]:


alpha = pow(1+(0.48+1.574*w - 0.176*w*w)*(1 - pow(Tr,0.5)),2)


# In[8]:


ac = 0.42728*pow(R*Tci, 2)/(Pci*pow(10,6))


# In[9]:


ai = ac*alpha


# In[11]:


k_data_file = open(r"k_data.csv","r")
csv_reader = csv.reader(k_data_file)

k_data =[]

for i in csv_reader:
    k_data.append(i[1:])

k_data = k_data[1:]


# In[13]:


a = 0
for i in range(1,12):
    for j in range(i,12):
        tmp = yi[i] * yi[j] * (( ai[i] * ai[j] )**0.5) * (1 - float(k_data[i][j]))
        a += tmp


# In[14]:


b = sum(yi*bi)


# In[16]:


P = R*T/(Vm - b) - a/(Vm*(Vm+b) + b*(Vm - b))


# In[17]:


print(P)


# In[ ]:




