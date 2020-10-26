# -*- coding: utf-8 -*-
"""
Created on %17/09/20

@author: %drylab iGEM UNIL (Pablo, Johanna, Julien)

"""
#Model of a killswitch type 2, whrere the action of the antitoxin occcurs at the protein level
#Check the modelling section of our wiki to understand the derivation of the equations
#Assuming that both repressors (xi, and p) are in excess and the amount in the cell will be the external amount 
#Equations from demonstration of hill equation igem tutorial week 2a 2/3
#Transcription is fast enough as compared to translation, so it was assumed to be at quasi-steady state.
#We had a call with A.vignoni and he suggested to introduce an equation on how the complex T-AT is evolving, and also another one for how the population is growing
#For the toxin they suggest this reference Synchronized cycles of bacterial lysis for in vivo delivery M. Omar Din, Tal Danino, Arthur Prindle, Matt Skalak, Jangir Selimkhanov, Kaitlin Allen, Ellixis Julio, Eta Atolia, Lev S. Tsimring, Sangeeta N. Bhatia & Jeff Hasty



#IMPORT PACKAGES
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

#DEFINE THE EQUATIONS
def killswitch_typeII(z, t):
    
    T     = z[0] #TOXIN
    AT    = z[1] #ANTITOXIN
    C     = z[2] #COMPLEX T-AT
    N     = z[3] #POPULATION
   
    ##################################  TOXIN EQUATION ######################################
    k1_T    = 0.78 #/min   #transcription rate of the toxin
    k2_T    = 1.78 #au (salis predictor/10000)     #translation rate of the toxin
    d1_T    = 0.231   #/min  #degradation rate of the toxin mRNA
    d2_T    = 0.0003   #/min #degradation rate of the toxin protein
    cn_T    = 2  #copy numbers of the toxin gene
    gamma_f = 1.67 #/(M*min)       #formation of the complex p_T + p_AT
    gamma_d = 0.029 #/(M*min)   #dissociation of the complex p_T + p_AT
    mu      = 0.085   #/min  #growth rate
    

    dprot_Tdt = k2_T*k1_T*cn_T/d1_T - d2_T*T - mu*T - gamma_f*T*AT + gamma_d*C
    
    
    
    ################################  ANTITOXIN EQUATION ####################################
    k1_AT  = 0.36  #/min #transcription rate of the antitoxin
    k2_AT  = 0.73  #translation rate of the antitoxin
    phosph = 0#repressor of the antitoxin expression 
    d1_AT  = 0.21      #degradation rate of the antitoxin mRNA 
    d2_AT  = 0.0003     #degradation rate of the antitoxin protein
    cn_AT  = 20 #copy numbers of the antitoxin gene # 2*5 copies of the plasmid
    kd_AT  = 1000 #molecules     #I don´t know this value 
    n_P    = 1 
    dprot_ATdt = k2_AT*k1_AT*cn_AT*kd_AT/(d1_AT*(
        kd_AT + phosph^n_P)) - d2_AT*AT - mu*AT - gamma_f*T*AT + gamma_d*C


    #################################  COMPLEX EQUATION #####################################
    d_Comp = 0     #also have to do research for this
    
    dCompdt = gamma_f*T*AT - gamma_d*C - d_Comp*C -mu*C

    
    
    ###############################  POPULATION DYNAMICS #####################################

    Nmax = 1000
    kmax = 10
    k0   = 2
    n_T=2
    dNdt = mu*N*(1 - N/Nmax) - N*(kmax*T**n_T)/(k0 + T**n_T) 
    
    #the first term of this equation is the logistic growth
    
    return[dprot_Tdt, dprot_ATdt, dCompdt, dNdt]


#initial conditions
z0 = [0,0,0,100]

#timepoints
t = np.linspace(0,200,500)

#INTEGRATE THE EQUATIONS IN LITTLE STEPS
z=odeint(killswitch_typeII, z0, t)

#PLOT
plt.plot(t,z[:,3],label = "population", color="purple", linewidth=3)
plt.plot(t,z[:,0],label = "toxin", color= "coral", linewidth=3)
plt.plot(t,z[:,1],label = "antitoxin", color="green",linewidth=3)
plt.xlabel("Time (min)", fontsize=14)
plt.legend(fontsize=14)
plt.ylim([0,1000])
plt.savefig("no_p_2_plasmid", dpi=500)
plt.show()














  
