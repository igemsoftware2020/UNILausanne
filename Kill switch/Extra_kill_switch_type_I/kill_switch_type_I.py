# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 8:30:31 2020

@author: pablo
"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
def repressilator(z, t):

#MODEL FOR KILLSIWTCH TYPE I: THE EFFECT OF THE ANTITOXIN IS AT THE MRNA LEVEL
#assuming that both repressors (xi, and p) are in excess and the amount in the cell will be the external amount
#also assuming that T-AT bind immediatelly
#equations from demonstration of hill equation igem tutorial week 2a 2/3
#also we assume that there is no basal expression
#hill coefficient would be 1 in this case, so we do not introduce n in the equation for now (but recheck this)
#we had a call with A.vignoni and he suggested to introduce an equation on how the complex T-AT is evolving, and also another one for how the population is growing
    
    m_T = z[0] #mRNA toxin
    p_T = z[1] #protein toxin
    m_AT = z[2] # mRNA antitoxin
    C = z[3] #complex mRNAtoxin-mRNAantitoxin
    N = z[4] #population
   
    ##### (If finally we use those values round them to 2 decimals)
    k1_1 = 2.157 #transcription rate of the toxin
    k2_1 = 0.826      #traduction rate of the toxin
    xi=0         #repressor of the toxin expression (we finally might not have it, so 0 for now)
    d1_1= 0.231      #degradation rate of the mRNA of the toxin
    d2_1= 0.0003      #degradation rate of the protein of the toxin itself
    cn_1= 3       #copy number of the toxin, I have 3 because we plan to have it repeated
    kd_1= 200     #disociation of the inhibitor of the inhibitor of the repressor for the toxin (as I said above it does not matter right now)
    gamma = 5      #formation of the complex m_T+m_AT
    gamma_m = 0.5  #dissociation of the complex m_T+m_AT
    mu = 0.085       #growth rate
    n=1
    
    dm_Tdt= k1_1*cn_1*kd_1/(kd_1+xi^n)-d1_1*m_T-gamma*m_T*m_AT-gamma_m*C-mu*m_T
    
    #
    
    dp_Tdt= k2_1*m_T - d2_1*p_T - mu*p_T
    
    #####
    k1_2 = 0.32
    p = 1000
    d1_2 = 0.231
    cn_2 = 20 #it will be ~20 for the repressilator and ~5 for the sponge
    kd_2 = 200 #I don´t know this value yet
    
    dm_ATdt=k1_2*cn_2*kd_2/(kd_2+p^n)-d1_2*m_AT-gamma*m_T*m_AT+gamma_m*C-mu*m_T


    ##### how the complex m_T-m_AT is working
    d_C = 0.5 #also have to do research for this
    
    dCdt=gamma*m_T*m_AT-gamma_m*C - d_C*C -mu*C

    #####
    Nmax = 10000
    kmax = 10
    k0 = 2
    dNdt = mu*N*(1-N/Nmax)-(kmax*p_T)/(k0*p_T)
    
    return[dm_Tdt, dp_Tdt, dm_ATdt, dCdt, dNdt]

  
#initial condition
z0=[20,20,20,0,200]#I made this up
#timepoints
t=np.linspace(0,200,500)


z=odeint(repressilator, z0, t)

plt.plot(t,z[:,0],label="mRNA_toxin")
plt.plot(t,z[:,1],label="p_toxin")
plt.plot(t,z[:,2],label="mRNA_antitoxin")
plt.plot(t,z[:,3],label="complex")
plt.plot(t,z[:,4],label="population")
plt.legend()
plt.xlabel("time")
plt.ylabel("protein expression")
plt.show()



