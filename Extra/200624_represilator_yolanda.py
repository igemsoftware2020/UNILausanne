# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 23:30:31 2020

@author: pablo
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
def repressilator(z, t):
    mRNA_lacI=z[0]
    p_lacI=z[1]
    mRNA_tetR=z[2]
    p_tetR=z[3]
    mRNA_cI=z[4]
    p_cI=z[5]
    #
    a_lacI,k_lacI,a_0,gamma,n=(216,1,0,1,5)
    dmRNA_lacIdt=a_lacI*(k_lacI**n/(k_lacI**n+p_cI**n)) + a_0 - gamma*mRNA_lacI
    #
    beta,delta=(5,5)
    dp_lacIdt=beta*mRNA_lacI-delta*p_lacI
    #
    a_tetR,k_tetR=(216,1)
    dmRNA_tetRdt=a_tetR*(k_tetR**n/(k_tetR**n+p_lacI**n)) + a_0 - gamma*mRNA_tetR
    #
    dp_tetRdt=beta*mRNA_tetR-delta*p_tetR
    #
    a_cI,k_cI=(216,1)
    dmRNA_cIdt=a_cI*(k_cI**n/(k_cI**n+p_tetR**n)) + a_0 - gamma*mRNA_cI
    #
    dp_cIdt=beta*mRNA_cI-delta*p_cI
    #
    return[dmRNA_lacIdt,dp_lacIdt,dmRNA_tetRdt,dp_tetRdt,dmRNA_cIdt,dp_cIdt]

  
#initial condition
z0=[0,0,0.1,0,0.2,0]#I made this up
#timepoints
t=np.linspace(0,50,500)


z=odeint(repressilator, z0, t)

plt.plot(t,z[:,1],label="p_lacI")
plt.plot(t,z[:,3],label="p_tetR")
plt.plot(t,z[:,5],label="p_cI")
plt.legend()
plt.xlabel("time")
plt.ylabel("protein expression")
plt.show()



