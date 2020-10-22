# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:02:13 2020

@author: pablo
"""


#I import the packages we will use
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
from matplotlib import colors
from statistics import stdev 
from statistics import mean


#################
# FUNCTIONS #
#################

#REPRESSILATOR
IPTG_0=0
def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    #
    n, nIPTG, Kd=( 2.4, 1, 10^-10)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n)
    #
    dp_cIdt = alpha/(1+p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]

        
#TO INITIATE OUR LIST
def initiation(number):
    for i in range(number):
        grid.append(initial_cell.copy())



def division(probability):
    #probability=np.random.normal(probablity, 0.1, 1)[0] #Uncoment this if you want to sample the probability from a normal distribution
    if random.uniform(0,1) <= probability: 
        new_rep=[x/2 for x in grid[cell].copy()]
        variability= 1
        grid[cell]=[x*variability for x in new_rep]
        cell_2=grid[cell].copy()
        cell_2=[x*(2-variability) for x in new_rep]
        grid.append(cell_2)
    
#How to get a global signal form all the cells    
def global_signal():
    p=[sum(x)*40/len(grid) for x in zip(*grid)]#we normalize dividing by grid len equivalent to dividing by od; *40 to plot the real amount of proteins (normalization in Elowitz et Leiber 2000)
    return(p[0],p[1],p[2])



# WE DEFINE THE PARAMETERS
grid=[]
time_limit=100;time=0
max_cells=1000 #this is our OD=0.4
cells_after_dilution=200 #this is our OD=0.2
div_prob=1
initial_cell=[0,4,0]#this aims to be the state of the three cells of the repressilator
# We set all the parameters through we want to loop:
alpha_= [0,50,100,150,200,250,300,350,400,450,500]
integration=[0.001,0.01,0.1,1,10,100]
# We create a thata frame to store teh results:
df = pd.DataFrame(0, columns = integration, index=alpha_)

#THERE WE GO!      
for alpha in alpha_:
    for integr in integration:
        time_limit=100;time=0
        initiation(1)
        proteins=[[],[],[]]
        while time<time_limit:
            for cell in range(len(grid)):
                  wi = odeint(repressilator, grid[cell], [0,integr] )
                  grid[cell]=wi[1]
                  division(div_prob)
            wa=global_signal()
            proteins[0].append(wa[0])
            proteins[1].append(wa[1])
            proteins[2].append(wa[2])
            if len(grid)>= max_cells:
                grid=random.sample(grid,cells_after_dilution)
                #when dividing cells should choose a random position, so this should still be a random sample of the population
            time+=1
            
        # If you want to plot the oscillations for each of the conditions in the loop
        plt.plot(range(time_limit),proteins[0],label="lacI")
        plt.plot(range(time_limit),proteins[1],label="tetR")
        plt.plot(range(time_limit),proteins[2],label="cI")
        plt.legend()
        plt.title("alpha: " + format(alpha) + " ; time_integration: " + format(integr))
        plt.xlabel("time (generations)")
        plt.ylabel("protein expression")
        #plt.savefig("data_2_4/fig_" + format(alpha) + format(integr) + ".png")
        plt.show()
        #to get an approximation of the amplitude towards the end, I just do the difference between the max and the min value of protein expression in the last 20 generations
        prot_1=max(proteins[0][80:])-min(proteins[0][80:])
        prot_2=max(proteins[1][80:])-min(proteins[1][80:])
        prot_3=max(proteins[2][80:])-min(proteins[2][80:])
        #We store teh value of the amplitude
        df[integr][alpha]=max(prot_1,prot_2,prot_3)
        #If we want to print in the terminal to keep tral
        print(alpha, "/", integr, "/",max(prot_1,prot_2,prot_3))

#######################################
# PLOTTING AND SAVING THE INFORMATION #
#######################################


df.to_csv("data_2_4/201018_alpha_integr.csv") #store wherever you want
data = pd.read_csv("data_2_4/201018_alpha_integr.csv") #load what you stored


#FINAL PLOT WE REPRESENT IN THE FIGURE IN THE WIKI
#instead of data.iloc[:,1:] you can just put df if you just run it
plt.contourf(data.iloc[:,1:], levels=[0,20,200,2000,20000,200000,2000000],colors=["black","beige","khaki","darkkhaki","olive","darkolivegreen"])
colores=["black","beige","khaki","darkkhaki","olive","darkolivegreen"]
plt.plot([0, 5], [4.32,4.32], ':', lw=3,color="white")
plt.plot([2.7,2.7], [0,10], ':', lw=3,color="white")
proxy = [plt.Rectangle((0,0),1,1,fc = pc,ec="black") for pc in colores]
plt.legend(proxy, ["range(0-2e1)", "range(2e1-2e2)", "range(2e2-2e3)", "range (2e3-2e4)", "range(2e4-2e5)","range(2e5-2e6)"],framealpha=0.5, loc='upper right')
plt.xticks(range(6),["1e-3","1e-2","1e-1","1e0","1e1","1e2","1e3"], rotation=0)
plt.yticks(range(11),alpha_)
plt.xlabel("integration_time (arbitrary units)")
plt.ylabel("alpha")
plt.title("Final amplitude of the oscillations (proteins/cell)")
plt.savefig('demo.png', transparent=True)
