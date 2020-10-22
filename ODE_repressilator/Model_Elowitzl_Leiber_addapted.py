# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 23:30:31 2020
Modified on Wed Jul 15 11:33:30 2020

@author: pablo
@made pretty: Johanna


"""
# COMENTS

# Here we used the model from Elowitz and Leiber 2000 and adapted it by introducign IPTG and aTc into the equations
# Note all the rescaling done in the model, we will have to rescale back protein and time before ploting
# IPTG "kidnaps" LacI and aTc "kidnaps" TetR
# We have introduced those compounds using the Hill Function thans to the tutorials of Imperial College iGEM TEAM 2020
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2715899/ #IPTG PAPER
# http://2016.igem.org/Team:William_and_Mary/Model #aTc IFO
# Nice quote:
# "[Mathematical models are] accurate descriptions of our pathetic thinking about nature" - James Black


# REQUIRED PACKAGES
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from statistics import mean



# PARAMETERS (those without units are adimentional)
alpha0 = 0.216
alpha  = 216
n      = 2
beta   = 0.2
IPTG_0 = 0
aTc_0  = 0
Kdi    = 1e-10 # M
Kda    = 3.6e-10 # M
ni     = 1
na     = 1


# We define the repressilartor as a function
def repressilator(z, t):
   
    mRNA_lacI = z[0]
    p_lacI    = z[1]
    mRNA_tetR = z[2]
    p_tetR    = z[3]
    mRNA_cI   = z[4]
    p_cI      = z[5]


    #LacI 

    dmRNA_lacIdt = alpha/(1 + p_cI**n) + alpha0 - mRNA_lacI
    dp_lacIdt    = - beta* (p_lacI - mRNA_lacI ) * (1- IPTG_0**ni/(Kdi**ni + IPTG_0**ni))
    
    #TetR 

    dmRNA_tetRdt = alpha/(1 + p_lacI**n) + alpha0 - mRNA_tetR
    dp_tetRdt    = - beta* (p_tetR- mRNA_tetR) * (1- aTc_0**na/(Kda**na + aTc_0**na))

    #CI
  
    dmRNA_cIdt   = alpha/(1+(p_tetR)**n)+ alpha0 - mRNA_cI
    dp_cIdt      = - beta* (p_cI -mRNA_cI)
    
    
    return[dmRNA_lacIdt, dp_lacIdt, dmRNA_tetRdt, dp_tetRdt, dmRNA_cIdt, 
           dp_cIdt]

  

#Initial conditions 
z0 = [0, 0, 1, 2, 0, 0] #mRNA, Protein: LacI, TetR, cI


#Timepoints
# Remember time is rescaled by ln(2)/2 min, so we undo the rescaling to polot
t  = np.linspace(0,500*(np.log(2)/2), 200) # we integrate 200 little stpept until 

t_to_plot=np.linspace(0,500*(np.log(2)/2), 6)
real_time=np.linspace(0, 500 ,6)
real_t=[int(round(x)) for x in real_time]
#Ploting

# We solve the equations and store their value.
z  = odeint(repressilator, z0, t)

# We create a list of proteins to plot
average_prot = []
for i in range(len(z)):
    average_prot.append(mean([z[i,1]*40, z[i,3]*40, z[i,5]]*40)) #times 40 because proteins are rescaled in units of KM

#We make the plot!
plt.plot(t,z[:,1]*40,label = "free LacI", color="tomato", lw=3)
plt.plot(t,z[:,3]*40,label = "free TetR", color="cornflowerblue", lw=3)
plt.plot(t,z[:,5]*40,label = "cI", color="y", lw=3)
plt.legend(fontsize=14, loc="upper right")
plt.xticks(t_to_plot,real_t, rotation=0, fontsize=12)
plt.ylim(0,3000)
plt.xlim(0,500*(np.log(2)/2))
plt.xlabel("Time (min)", fontsize=16)
plt.ylabel("Molecules per cell", fontsize=16)
plt.show()

# In case you want to create a plot per step and then merge them in a video or gif, uncoment and run the following
"""
for i in range(200):

    plt.plot(t[:i],z[:,1][:i]*40,label = "LacI", color="tomato", lw=3)
    plt.plot(t[:i],z[:,3][:i]*40,label = "TetR", color="cornflowerblue", lw=3)
    plt.plot(t[:i],z[:,5][:i]*40,label = "cI", color="blue", lw=3)
    #plt.plot(t, average_prot, label = "Average protein")
    plt.legend(fontsize=14, loc="upper right")
    plt.xticks(np.linspace(0, 34.657*5,6),real_t, rotation=0)
    plt.ylim(0,2500)
    plt.xlim(0,34.657*5)
    plt.xlabel("Time (min)", fontsize=14)
    plt.ylabel("Molecules per cell", fontsize=14)
    plt.savefig("WHEREVER/fig_" + format(i) +".png", transparent=True)
    plt.show()

#plt.plot(z[:,0], z[:,1])
    """