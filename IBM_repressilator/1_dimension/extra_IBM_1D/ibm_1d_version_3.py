# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 10:18:59 2020

@author: pablo
"""
#this aims to be a 1d simplification of the individual based model to assess how noise during bacterial division can affect the repressilator.
#ok, so en the previous version I realised that if the noise I introduced was divide or not, the more noise the fewer the speed of division and thus the speed of oscillation. So now I will try that teh bias goes in both sides (a cell can divide, not divide or divide twice).



#I import the packages I will use
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors
from statistics import stdev 
from statistics import mean

#noise=0.1
#s = np.random.choice([0,1,2],50,p=[noise/2,(1-noise),noise/2])
#count, bins, ignored = plt.hist(s, 3, normed=True)
#plt.show()



#################
# Repressilator #
#################
IPTG_0=0
def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    #
    alpha,n, nIPTG, Kd=(100,2.4,2, 10^-10)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n)
    #
    dp_cIdt = alpha/(1+p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]




##########
# IBM_1D #
##########
#noise_div=0.1 #from 0 to 1 (this should be somehow the inverse to my division probability in my 2D model)
prob_div=1 #somehow noise will be 1-prob_div
initial_cell={"div":0,"rep":[0,4,0]}#this aims to be the state of the three cells of the repressilator

doubling_time=35 #we assume that the doubling time of bacteria is this in minutes
generations=50
time_limit=generations*doubling_time;time=0
max_cells=1000 #this is our OD=0.4
cells_after_dilution=200 #this is our OD=0.2


grid=[]

def initiation(number):
    for i in range(number):
        grid.append(initial_cell.copy())
    

def division(prob):
    if random.uniform(0,1) <= prob:
        new_rep=[x/2 for x in grid[cell]["rep"].copy()]
        variability=1# np.random.normal(1, 0.2, 1)[0]
        grid[cell]["rep"]=[x*variability for x in new_rep]
        grid[cell]["div"]-=1
        cell_2=grid[cell].copy()
        cell_2["rep"]=[x*(2-variability) for x in new_rep]
        grid.append(cell_2)
    
def global_signal():
    p=[sum(x)/len(grid) for x in zip(*grid)]#we normalize dividing by grid len equivalent to dividing by od
    return(p[0],p[1],p[2])
    #maybe I have to divide this by the number of cells

def all_values():
    p=[[],[],[]]
    for i in grid:
        p[0].append(i["rep"][0])
        p[1].append(i["rep"][1])
        p[2].append(i["rep"][2])
    return(p)


#NOW LET´S GO WITH THE LOOP =^_^= 
initiation(2) #I have to start at least with 2 cells to calculate the sdanderd deviation
proteins=[[],[],[]]
s_dev=[[],[],[]]
while time<time_limit:
    for cell in range(len(grid)):
        grid[cell]["div"]+=1
        wi = odeint(repressilator, grid[cell]["rep"], [0,0.5] )
        grid[cell]["rep"]=wi[1]
        division(prob_div)
    wa=all_values()
    proteins[0].append(mean(wa[0]))
    proteins[1].append(mean(wa[1]))
    proteins[2].append(mean(wa[2]))
    s_dev[0].append(stdev(wa[0]))
    s_dev[1].append(stdev(wa[1]))
    s_dev[2].append(stdev(wa[2]))
    if len(grid)>= max_cells:
        grid=random.sample(grid,cells_after_dilution)
        #when dividing cells should choose a random position, so this should still be a random sample of the population
    time+=doubling_time
    
      
plt.plot(np.linspace(0,time_limit,generations),proteins[0],label="lacI", color="tomato")
plt.plot(np.linspace(0,time_limit,generations),proteins[1],label="tetR", color="cornflowerblue")
plt.plot(np.linspace(0,time_limit,generations),proteins[2],label="cI", color="y")
plt.errorbar(np.linspace(0,time_limit,generations),proteins[0],s_dev[0], lw=0.5, color="tomato")
plt.errorbar(np.linspace(0,time_limit,generations),proteins[1],s_dev[1], lw=0.5, color="cornflowerblue")
plt.errorbar(np.linspace(0,time_limit,generations),proteins[2],s_dev[2], lw=0.5, color="y")
plt.legend()
plt.xlabel("Time (min)")
plt.ylabel("Average molecules per cell ± SD")
plt.show()

"""
plt.plot(np.linspace(0,time_limit,generations),proteins[0],label="lacI", color="tomato")
plt.plot(np.linspace(0,time_limit,generations),s_dev[0],label="sd_lacI",color="cornflowerblue")
plt.xlabel("time (min)")
plt.legend()
plt.show()

division_track=[]
for cell in range(len(grid)):
    division_track.append(grid[cell]["div"])
values, counts = np.unique(division_track, return_counts=True)
counts_freq=[x/len(division_track) for x in counts]
plt.vlines(values, 0, counts_freq, color='C0', lw=20)
plt.ylabel("frequency")
plt.show()
# optionally set y-axis up nicely
#plt.ylim(0, max(counts_freq) * 1
"""