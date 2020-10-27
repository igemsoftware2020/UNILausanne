#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 22:43:16 2020

@author: pablo
"""

#MODEL TO INVESTIGATE HOW THE EQUATIONS IN OUR ODE MODEL LOOK WITHOUT CELL DIVISION

#I import the packages I will use
import numpy as np
import random
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import odeint
from matplotlib import colors




#################
# Repressilator #
#################


def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    IPTG_0=0.00001
    alpha,n, nIPTG, Kd=(216,2.4,1, 1.4e-6)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd + IPTG_0**nIPTG))#-0.0001*p_lacI
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n) #- 0.0001*p_tetR
    #
    dp_cIdt=alpha/(1+p_tetR**n) #- 0.0001* p_cI
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


##########################
# Individual Based Model #
##########################
IPTG_0=0
death_rate=0
grid_size=100 #here I set the grid size, which will be the number of rows and columsn of a square grid
time_steps=3091; time_step=0#numer of rounds and the initial round for the while loop

cell_0={"div":0, "rep":[10,1,1]} 

protein1=[]
protein2=[]
protein3=[]

while time_step<time_steps: #the simulation will run until the number of rounds previously set
     wi = odeint(repressilator, cell_0["rep"], [0,50] )
     cell_0["rep"]=wi[1]
     protein1.append(wi[1][0])
     protein2.append(wi[1][1]) 
     protein3.append(wi[1][2])            
                        
     time_step+=1

time=range(time_steps)
plt.plot(time, protein1, label = "p_LacI")
plt.plot(time, protein2, label = "p_TetR")
plt.plot(time, protein3, label = "p_cI")
#plt.plot(time, media, label = "Average protein")
plt.legend()
plt.xlabel("Time")
plt.ylabel("Protein expression")
plt.show()