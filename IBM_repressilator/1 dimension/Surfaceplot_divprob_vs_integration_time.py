# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 12:02:13 2020

@author: pablo
"""

#We did not have time to fully clean this script, still it works, for a better organization check "Surgaceplot_alpha_vs_integration_time.py" as it is basically the same and it is organised in a better way there and has clearer comments.

#I import the packages I will use
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.integrate import odeint
from matplotlib import colors
from statistics import stdev 
from statistics import mean

#noise=0.1
#s = np.random.choice([0,1,2],50,p=[noise/2,(1-noise),noise/2])
#count, bins, ignored = plt.hist(s, 3, normed=True)
#plt.show()

prob_division=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
integration=[0.001,0.01,0.1,1,10,100]
df = pd.DataFrame(0, columns = integration, index=prob_division)

for div_prob in prob_division:
    for integr in integration:
        #################
        # Repressilator #
        #################
        IPTG_0=0
        def repressilator(z, t):
            p_lacI=z[0]
            p_tetR=z[1]
            p_cI=z[2]
            #
            alpha,n, nIPTG, Kd=( 216, 2.4, 1, 10^-10)
            dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
            
            #
            dp_tetRdt = alpha/(1+p_lacI**n)
            #
            dp_cIdt = alpha/(1+p_tetR**n)
            return[dp_lacIdt,dp_tetRdt,dp_cIdt]
        
        
        
        
        ##########
        # IBM_1D #
        ##########
        #div_prob=1 #from 0 to 1 (this should be somehow the inverse to my division probability in my 2D model)
        
        initial_cell=[0,10,0]#this aims to be the state of the three cells of the repressilator
        
        time_limit=100;time=0
        max_cells=1000 #this is our OD=0.4
        cells_after_dilution=200 #this is our OD=0.2
        
        
        grid=[]
        
        def initiation(number):
            for i in range(number):
                grid.append(initial_cell.copy())
        
        
        
        def division(wiii):
            prob=np.random.normal(wiii, 0.1, 1)[0] #I want to try what sara suggested, lets keep sd=0.1 for now
            if random.uniform(0,1) <= prob:
                new_rep=[x/2 for x in grid[cell].copy()]
                variability= 1# np.random.normal(1, 0.2, 1)[0] #or also from a beta distribution as we did later (between 0 and 1)
                grid[cell]=[x*variability for x in new_rep]
                cell_2=grid[cell].copy()
                cell_2=[x*(2-variability) for x in new_rep]
                grid.append(cell_2)
            
                
        
        
        def global_signal():
            p=[sum(x)*40/len(grid) for x in zip(*grid)]#we normalize dividing by grid len equivalent to dividing by od; *40 to plot the real amount of proteins (normalization in elowitz)
            return(p[0],p[1],p[2])
            #maybe I have to divide this by the number of cells
        
        def all_values():
            p=[[],[],[]]
            for i in grid:
                p[0].append(i[0])
                p[1].append(i[1])
                p[2].append(i[2])
        
        
        #NOW LET´S GO WITH THE LOOP =^_^= 
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
            
              
        plt.plot(range(time_limit),proteins[0],label="lacI")
        plt.plot(range(time_limit),proteins[1],label="tetR")
        plt.plot(range(time_limit),proteins[2],label="cI")
        plt.legend()
        plt.title("div_prob: " + format(div_prob) + " ; time_integration: " + format(integr))
        plt.xlabel("time (generations)")
        plt.ylabel("protein expression")
        #plt.savefig("data_2_4/fig_" + format(div_prob) + format(integr) + ".png")
        plt.show()
        #to get an approximation of the amplitude towards the end, I just do the difference between the max and the min value of protein expression in the last 20 generations
        prot_1=max(proteins[0][80:])-min(proteins[0][80:])
        prot_2=max(proteins[1][80:])-min(proteins[1][80:])
        prot_3=max(proteins[2][80:])-min(proteins[2][80:])
        
        df[integr][div_prob]=max(prot_1,prot_2,prot_3)
        print(div_prob, "/", integr, "/",max(prot_1,prot_2,prot_3))

#At the end we did not use that, but just in case
sns.heatmap(df, annot=True, linewidths=.5, cmap= "tab10", center=1)
#in this case for example would be better to plot the log to show better
sns.heatmap(np.log(df), annot=True, linewidths=.5, center=0)
#http://seaborn.pydata.org/generated/seaborn.heatmap.html
df.to_csv("data_2_4/201019_alpha_div.csv")
# check plt.contourf to represent, also surface plots (you give the matrix and they interpolate to be smoth)
data = pd.read_csv("data_2_4/201019_alpha_div.csv") 
#plt.contourf(data.iloc[:,1:], levels=[0,1,10,100,1000,10000,1000000],colors=["black","papayawhip","moccasin","khaki","goldenrod","darkgoldenrod"])
#data = pd.read_csv("data/201014_alpha_integr.csv") 

plt.contourf(data.iloc[:,1:], levels=[0,20,200,2000,20000,200000,2000000],colors=["black","beige","khaki","darkkhaki","olive","darkolivegreen"])
colores=["black","beige","khaki","darkkhaki","olive","darkolivegreen"]
#plt.plot([0, 5], [4.32,4.32], ':', lw=3,color="white")
#plt.plot([2.7,2.7], [0,10], ':', lw=3,color="white")
proxy = [plt.Rectangle((0,0),1,1,fc = pc,ec="black") for pc in colores]
plt.legend(proxy, ["range(0-2e1)", "range(2e1-2e2)", "range(2e2-2e3)", "range (2e3-2e4)", "range(2e4-2e5)","range(2e5-2e6)"],framealpha=0.5, loc='upper right')
plt.xticks(range(6),["1e-3","1e-2","1e-1","1e0","1e1","1e2","1e3"], rotation=0)
plt.yticks(range(10),prob_division)
plt.xlabel("integration_time (arbitrary units)")
plt.ylabel("alpha")
plt.title("Final amplitude of the oscillations (proteins/cell)")
plt.savefig('new.png', transparent=True)
