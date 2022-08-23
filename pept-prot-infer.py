#!/usr/bin/env python
# coding: utf-8

# In[18]:


import numpy as np # better arrays than inbuilt arrays
import matplotlib.pyplot as plt # to plot stuff

import pandas as pd #for DataFrame tables
from IPython.display import display #to display dfs more nicely

import scipy.stats
from scipy.stats import norm
import statistics


# ## Generation of random distributions

# In[26]:


# Creation of random data points from multiple Gaussians ki
mean1, std1 = 120, 15 # mean and standard deviation
mean2, std2 = 170, 15

k1 = np.random.normal(mean1, std1, 15) #creates array with values created through Gaussian
k2 = np.random.normal(mean2, std2, 148)

x = np.concatenate([k1, k2])


# In[27]:


# Calculating p_xi_given_zj (this is what Matt is working on with the simulated data)

pdf_probability_k1 = scipy.stats.norm.pdf(x, loc=mean1, scale=std1)
pdf_probability_k2 = scipy.stats.norm.pdf(x, loc=mean2, scale=std2)

p_xi_given_zj = np.vstack((pdf_probability_k1,pdf_probability_k2))
display(pd.DataFrame(p_xi_given_zj))

print(np.argmax(p_xi_given_zj, axis=0)) # for each column of P(Xi|Zj), the most likely Peptide is returned


# In[28]:


# plotting histograms
nbins = 50
plt.hist(k1, label = "Peptide 0", bins=nbins, alpha=0.3, density=True, color="orange") # alpha=transparency, density=True normalises to 1 
plt.hist(k2, label = "Peptide 1", bins=nbins, alpha=0.3, density=True, color="green")

# PDF plot
xmin, xmax = plt.xlim() #finds lower and upper bounds of histogram data
x = np.linspace(start=xmin, stop=xmax, num=100) #num is the number of returned data points - the more points, the finer the fit is plotted
p1 = norm.pdf(x, mean1, std1)
p2 = norm.pdf(x, mean2, std2)
  
plt.plot(x, p1, linewidth=2, color = "orange", label = "Gauss function k1: mean = {:.2f}, STD = {:.2f}".format(mean1, std1))
plt.plot(x, p2, linewidth=2, color = "green", label = "Gauss function k2: mean = {:.2f}, STD = {:.2f}".format(mean2, std2))

plt.legend(loc='upper right')
plt.title("PDFs of dwarves and humans")

plt.show()


# ## First half
# First iteration:
# - assume P(Zi)=1/m for all i
# 
# Real data: a bunch of dye seq data points Xi.
# 
# Goal: Figure out probability of real peptide seqs Zj given all those measured Xi.
# ____
# However, I am starting with the list of P(Xi|Zj) that Matt will provide me, instead of a list of P(Xi). In fact, I never even use P(Xi), which is a little weird.
#     

# In[22]:


display(pd.DataFrame(p_xi_given_zj))

m = p_xi_given_zj.shape[0] # Number of rows = number of datapoints (dye seqs) x
n = p_xi_given_zj.shape[1] # Number of columns = number of peptides z
p_zj_initial = 1/n #initial approximation: all zj equally likely, to jumpstart first iteration

p_zj = np.full(n, p_zj_initial) 


# In[23]:


p_zj_given_xi = np.full((n, m), 0, dtype=float) #rows = j (peptides), columns = i (data)
denominator = 0
loopcounter = 0

while loopcounter <= 30:
    for i, row in enumerate(p_xi_given_zj): # Calculating/Updating P(Zj|Xi)
        # print("ROW:", i)
        for j, cell in enumerate(row):
            # print("COLUMN:", j)
            numerator = cell * p_zj[j]
            # print("numerator:", numerator)

            for l, cell in enumerate(p_zj):
                # print("cell", i, l, p_xi_given_zj[i, l], end="")
                # print(" * zl", p_zj[l])
                denominator = denominator + p_xi_given_zj[i, l] * p_zj[l]
            p_zj_given_xi[j][i] = numerator/denominator
            denominator = 0                                
    display(pd.DataFrame(p_zj_given_xi))
#    print(np.amax(p_zj_given_xi, axis=1)) #reports max value from each row
    print(np.argmax(p_zj_given_xi, axis=1)) # reports index of max value from each row

    for j, element in enumerate(p_zj): #updating the expectation value of Zi
        p_zj[j] = p_zj_given_xi[j].sum()/n
    display(pd.DataFrame(p_zj))
    
    loopcounter = loopcounter + 1


# ## Bootstrapping part
1. Random subtable
df = df.sample(frac=0.50,axis='columns', replace=true) #randomly selects fraction of columns, with replacement

2. EM, save estimate for each parameter
3. create a confidence interval each parameter by taking the middle x percent
# In[ ]:




