# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 15:42:38 2020

@author: pablo
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#We tried to apply the equations from boada et al. 2020 and the tutorials from Alejandro Vignoni to obtain the equations od the repressialtor
#all the parameters here are probably wrong, I just used those from Yolanda webpage to obtain oscillations
#we have to improve the parameters (check iGEM 2007 and 2011 zurich, and igem 2009 aberdeen).
#also try the RSB calculator of Salis lab to infer the translation rate


IPTG=0
def repressilator(z, t):
    lacI=z[0]
    tetR=z[1]
    cI=z[2]
    
    #
    k1_mlacI = 2
    k2_lacI = 5
    CN = 5
    d_mlacI = 0.231 #for now I assume is the same for all of them #this probably needs to be adapted to introduce dilution by cell growth
    alpha_pcI = 0 #let´s assume no leakage, if it is going to be finally like that we can remove it from the equations
    kd_pcI = 1 #zurich2007
    n = 5
    d_lacI = 1 #this probably needs to be adapted to introduce dilution by cell growth

    
    
    dlacIdt = (k1_mlacI*k2_lacI*CN/d_mlacI)*(alpha_pcI+(1-alpha_pcI)*kd_pcI**n/(kd_pcI**n + cI**n)) - d_lacI*lacI

    #
    k1_mtetR = 2
    k2_tetR = 5
    d_mtetR = 0.231
    alpha_placI = 0
    kd_placI = 1
    kd_lacI_IPTG =1.3e-6 #ethz 2007
    d_tetR = 1 #ethz 2007
    
    dtetRdt =(k1_mtetR*k2_tetR*CN/d_mtetR)*(alpha_placI+((1-alpha_placI)*((kd_placI*CN)**n)*(kd_lacI_IPTG+IPTG)**n)/(((kd_placI*CN)**n)*((kd_lacI_IPTG+IPTG)**n)+(kd_lacI_IPTG*lacI)**n)) - d_tetR*tetR
    
    #
    k1_mcI = 2
    k2_cI = 5
    d_mcI=0.231
    alpha_ptet=0    
    kd_ptet=1
    d_cI = 1 
    dcIdt = (k1_mcI*k2_cI*CN/d_mcI)*(alpha_ptet+(1-alpha_ptet)*kd_ptet**n/(kd_ptet**n + tetR**n)) -d_cI*cI
    
    return[dlacIdt,dtetRdt,dcIdt]


#initial condition
z0=[0.1,0,0]#I made this up
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