# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:51:31 2020

@author: pablo
"""
#this version of a minimal repressilator was taken from a professor in youtube (https://www.youtube.com/watch?v=xNNxlsY-F-s&t=3427s)
#it can be obtained from the model i the paper mathematical biosciences (page 8) if we adjust the units to remove the K as in the model of the repressilator of Elowitz et al 2000 and we do not consider the leakiness of the promoter. 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

IPTG_0=0
def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    
    #
    alpha,n, nIPTG, Kd=(216,3,2, 10^-10)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd + IPTG_0**nIPTG))-p_lacI
    #
    dp_tetRdt = alpha/(1+p_lacI**n)-p_tetR
    #
    dp_cIdt=alpha/(1+p_tetR**n)-p_cI
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


#initial condition
z0=[0,10,1]#I made this up
#timepoints
t=np.linspace(0,20,500)


z=odeint(repressilator, z0, t)

plt.plot(t,z[:,0],label="p_lacI", color="tomato", lw=3)
plt.plot(t,z[:,1],label="p_tetR", color="cornflowerblue",lw=3)
plt.plot(t,z[:,2],label="p_cI", color="y",lw=3)
plt.legend(fontsize=12, loc="upper right")
plt.xlabel("Time (arbitrary units)",fontsize=14)
plt.ylabel("Molecules per cell",fontsize=14)
plt.show()