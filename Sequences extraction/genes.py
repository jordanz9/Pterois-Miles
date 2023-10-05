# -*- coding: utf-8 -*-
'''
@author: jordanz

This python script generates a plot that maps the lionfish genes based on coordinates provided in the 'toxin_genes.txt' file.

'''

import matplotlib.pyplot as plt
import matplotlib.patches as patch
import numpy as np

plt.style.use('default')
fig = plt.figure(figsize = (16,9))
ax = fig.add_subplot(222)
plt.yticks([])


data = np.loadtxt("toxin_genes.txt", delimiter='\t', dtype = "str")

start =[]
end = []

for i in range(1,len(data)):
    start.append(float(data[i][1]))
    end.append(float(data[i][2]))
    
toprint = []

for i in range(len(start)):
    toprint.append(tuple([start[i],end[i]]))
    
    
minimum = start[0]
maximum = end[-1]

colors = ['maroon','lightcoral','orange','yellowgreen','g','c']
k = 0

for i in range(len(start)):
    
    if i % 3 == 0 and i != 0 :
        k+=1
        
    rect = patch.Rectangle((start[i],0), end[i] - start[i], height = 0.5,
                           color = colors[k])

    ax.add_patch(rect) 
    

    
    print(k)

plt.xlim([minimum - 5000, maximum + 5000])  
plt.ylim([-1,1.5])   
ax.set_title("Mapping Genes on Contig 192 Pilon")

plt.show()
