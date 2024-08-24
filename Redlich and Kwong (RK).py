#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import csv


# In[17]:


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


# In[6]:


alpha = pow(Tr,-0.5)


# In[7]:


ac = 0.42728*pow(R*Tci, 2)/(Pci*pow(10,6))


# In[8]:


ai = ac*alpha


# In[20]:





# In[16]:


b = sum(yi*bi)
print(b)


# In[15]:


a = pow(sum(yi*pow(ai,0.5)),2)
print(a)


# In[18]:


P = (R*T/(Vm - b)) - (a/(Vm*(Vm+b)))


# In[19]:


print(P)


# In[ ]:




