# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:18:59 2020

@author: pablo

This aims to be a 1d simplification of the individual based model to assess how noise during 
bacterial division can affect the repressilator.
"""

#One of the initial versions of our IBM 1D

import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
#from matplotlib import colors

#################
# Repressilator #
#################
IPTG_0=0
def repressilator(z, t):
    p_lacI = z[0]
    p_tetR = z[1]
    p_cI   = z[2]
    
    #
    alpha,n, nIPTG, Kd = (100, 3,2, 10^-10)
    dp_lacIdt = (alpha/(1 + p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1 + p_lacI**n)
    #
    dp_cIdt = alpha/(1 + p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


##########
# IBM_1D #
##########

initial_cell = [0,4,0]       #State of the three cells of the repressilator
prob_div     = 0.2
time_limit   = 200
time         = 0
max_cells    = 1000          #OD=0.4
cells_after_dilution = 200   #OD=0.2

grid=[]

def initiation(number):
    for i in range(number):
        grid.append(initial_cell.copy())

def division(prob):
    if random.uniform(0,1) <= prob_div:
        new_rep = [x/2 for x in grid[cell].copy()]
        variability = 1     #np.random.normal(1, 0.2, 1)[0]
        grid[cell] = [x*variability for x in new_rep]
        grid.append([x*(2-variability) for x in new_rep])


def global_signal():
    p = [sum(x)/len(grid) for x in zip(*grid)]
    #Normalize by dividing by grid_len, equivalent to dividing by OD
    
    return(p[0],p[1],p[2])
    #maybe I have to divide this by the number of cells

#NOW LET´S GO WITH THE LOOP =^_^= 
initiation(1)
proteins = [[],[],[]]
while time < time_limit:
    for cell in range(len(grid)):
          wi = odeint(repressilator, grid[cell], [0,0.5] )
          grid[cell] = wi[1]
          division(prob_div)
    
    wa = global_signal()
    proteins[0].append(wa[0])
    proteins[1].append(wa[1])
    proteins[2].append(wa[2])
    
    if len(grid) >= max_cells:
        grid = random.sample(grid, cells_after_dilution)
        #since we are takin a random sample should be ok
    
    time += 5
    
      
plt.plot(np.linspace(0,200,40),proteins[0],label = "lacI")
plt.plot(np.linspace(0,200,40),proteins[1],label = "tetR")
plt.plot(np.linspace(0,200,40),proteins[2],label = "cI")
plt.legend()
plt.xlabel("time")
plt.ylabel("protein expression")
plt.show()
