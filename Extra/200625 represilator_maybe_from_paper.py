# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 23:32:23 2020

@author: pablo
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    #
    gamma_lacI,k_lacI,n=(20,0.1,3)
    dp_lacIdt=gamma_lacI*(k_lacI**n/(k_lacI**n+p_cI**n))-p_lacI   
    #
    gamma_tetR,k_tetR=(20,0.1)
    dp_tetRdt=gamma_tetR*(k_tetR**n/(k_tetR**n+p_lacI**n))-p_tetR
    #
    gamma_cI,k_cI=(20,0.1)
    dp_cIdt=gamma_cI*(k_cI**n/(k_cI**n+p_tetR**n))-p_cI
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


#initial condition
z0=[0,1,1]#I made this up
#timepoints
t=np.linspace(0,50,500)


z=odeint(repressilator, z0, t)

plt.plot(t,z[:,0],label="p_lacI")
plt.plot(t,z[:,1],label="p_tetR")
plt.plot(t,z[:,2],label="p_cI")
plt.legend()
plt.xlabel("time")
plt.ylabel("protein expression")
plt.show()