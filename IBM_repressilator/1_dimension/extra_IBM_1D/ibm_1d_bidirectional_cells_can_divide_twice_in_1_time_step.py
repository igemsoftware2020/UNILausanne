# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:18:59 2020

@author: pablo



"""
#This aims to be a 1d simplification of the individual based model to assess how noise during 
#bacterial division can affect the repressilator. Ok, so en the previous version I realised that
#if the noise I introduced was divide or not, the more noise the fewer the speed of division and 
#thus the speed of oscillation. So now I will try that teh bias goes in both sides (a cell can 
#divide, not divide or divide twice).

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors


#noise=0.1
#s = np.random.choice([0,1,2],50,p=[noise/2,(1-noise),noise/2])
#count, bins, ignored = plt.hist(s, 3, normed=True)
#plt.show()



#################
# Repressilator #
#################
IPTG_0 = 0
def repressilator(z, t):
   
    p_lacI = z[0]
    p_tetR = z[1]
    p_cI   = z[2]
    
    #
    alpha,n, nIPTG, Kd = (100, 3, 2, 10^-10)
    dp_lacIdt = (alpha/(1 + p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1 + p_lacI**n)
    #
    dp_cIdt = alpha/(1 + p_tetR**n)
    
    return[dp_lacIdt, dp_tetRdt, dp_cIdt]




##########
# IBM_1D #
##########

noise_div    = 0  #Inverse to my division probability in my 2D model)
initial_cell = [0,2,0]      #Quantity proteins
time_limit   = 100
time         = 0
max_cells    = 1000         #OD=0.4
cells_after_dilution = 500  #OD=0.2 

grid=[]

def initiation(number):
    for i in range(number):
        grid.append(initial_cell.copy())

def division(times):
    
    if times == 1:
        new_rep     = [x/2 for x in grid[cell].copy()] 
        grid[cell]  = new_rep
        #variability= 1 #np.random.normal(1, 0.2, 1)[0] 
        grid.append(new_rep)

    
    elif times == 2:
        new_rep     = [x/3 for x in grid[cell].copy()] 
        grid[cell]  = new_rep
        grid.append(new_rep)
        grid.append(new_rep)

def global_signal():
    p = [sum(x)/len(grid) for x in zip(*grid)]
    #We normalize by dividing by the grid_len, equivalent to dividing by OD
    
    return(p[0], p[1], p[2])

#NOW LET´S GO WITH THE LOOP =^_^= 
    
initiation(1)
proteins = [[],[],[]]
while time < time_limit:
    for cell in range(len(grid)):
          wi         = odeint(repressilator, grid[cell], [0,0.2])
          grid[cell] = wi[1]
          s        = np.random.choice([0,1,2],1,p = [noise_div/2,(1 - noise_div),noise_div/2])
          division(s[0])
    
    wa = global_signal()
    proteins[0].append(wa[0])
    proteins[1].append(wa[1])
    proteins[2].append(wa[2])
    
    if len(grid) >= max_cells:
        grid = random.sample(grid,cells_after_dilution)
        #when dividing cells should choose a random position, so this should still be a random 
        #sample of the population
    
    time += 1
    
      
plt.plot(range(time_limit),proteins[0],label = "lacI")
plt.plot(range(time_limit),proteins[1],label = "tetR")
plt.plot(range(time_limit),proteins[2],label = "cI")
plt.legend()
plt.xlabel("time")
plt.ylabel("protein expression")
plt.show()








