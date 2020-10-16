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
import scipy
from scipy.stats import shapiro 
from scipy.integrate import odeint
from matplotlib import colors
from statistics import stdev 
from statistics import mean
from statistics import variance
import copy

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
    alpha,n, nIPTG, Kd=(216,3,2, 10^-10)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n)
    #
    dp_cIdt = alpha/(1+p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]




##########
# IBM_1D #
##########

##Here I define the parameters
prob_div=1 #somehow noise will be 1-prob_div
initial_cell={"div":0,"rep":[0,4,0]}#this aims to be the state of the three cells of the repressilator
time_limit=100; time=0 #to plot autocorrelation I use 200 time steps instead
max_cells=1000 #this is our OD=0.4
cells_after_dilution=200 #this is our OD=0.2
cells_to_plot=10

grid=[]


##Funcions I need to use
def initiation(number):
    for i in range(number):
        grid.append(initial_cell.copy())
    

def division(prob):
    if random.uniform(0,1) <= prob: #here you can put "np.random.normal(prob,0.1):" instead of "prob:"
        new_rep=grid[cell]["rep"].copy()
        variability= 0.5 #np.random.beta(50,50)#I choose this instead of a normal distribution because a normal distribution potentially goes from minus infinity to infinity. If we don´t want variability we just put 0.5
        grid[cell]["rep"]=[x*variability for x in new_rep]
        grid[cell]["div"]+=1
        cell_2=grid[cell].copy()
        cell_2["rep"]=[x*(1-variability) for x in new_rep]
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
initiation(cells_to_plot) #I have to start at least with 2 cells to calculate the sdanderd deviation
proteins=[[],[],[]]#value of the three proteins
sem=[[],[],[]]#standard error of the mean
sd=[[],[],[]]#standard deviation
dict_to_plot={}
for i in range(cells_to_plot):
    dict_to_plot["lacI{}".format(i)]=[]
    dict_to_plot["tetR{}".format(i)]=[]
    dict_to_plot["cI{}".format(i)]=[]
while time<time_limit:
    grid_copy=copy.copy(grid)
    for cell in range(len(grid_copy)):
        wi = odeint(repressilator, grid[cell]["rep"], [0,5] )
        grid[cell]["rep"]=wi[1]
        division(prob_div)
    
    #we store all the information we need
    wa=all_values() #this would be a list with three list each of them has tha values for the amount of one of the proteins in all the cells
    proteins[0].append(mean(wa[0]))
    proteins[1].append(mean(wa[1]))
    proteins[2].append(mean(wa[2]))
    sem[0].append(scipy.stats.sem(wa[0]))
    sem[1].append(scipy.stats.sem(wa[1]))
    sem[2].append(scipy.stats.sem(wa[2]))
    sd[0].append(stdev(wa[0]))
    sd[1].append(stdev(wa[1]))
    sd[2].append(stdev(wa[2]))
                  
    for i in range(cells_to_plot):
        dict_to_plot["lacI{}".format(i)].append(grid[i]["rep"][0])
        dict_to_plot["tetR{}".format(i)].append(grid[i]["rep"][1])
        dict_to_plot["cI{}".format(i)].append(grid[i]["rep"][2])
    if len(grid)>= max_cells:
        grid=grid[:cells_after_dilution]# Usually I would take a random sample, but here I cannot as I want to keep track of the same cells in every time step
    time+=1

#PLOT +- SEM
time_to_plot=np.linspace(0,time_limit,time_limit)
plt.plot(time_to_plot,proteins[0],label="lacI", color="tomato")
plt.fill_between(time_to_plot,np.asarray(proteins[0])-np.asarray(sem[0]) , np.asarray(proteins[0])+ np.asarray(sem[0]), color="tomato", alpha=0.2, edgecolor='none')
plt.plot(time_to_plot,proteins[1],label="tetR", color="c")
plt.fill_between(time_to_plot,np.asarray(proteins[1])-np.asarray(sem[1]) , np.asarray(proteins[1])+ np.asarray(sem[1]), color="c", alpha=0.2, edgecolor='none')
plt.plot(time_to_plot,proteins[2],label="cI", color="y")
plt.fill_between(time_to_plot,np.asarray(proteins[2])-np.asarray(sem[2]) , np.asarray(proteins[2])+ np.asarray(sem[2]), color="y", alpha=0.2, edgecolor='none')
plt.legend(fontsize=12,framealpha=0.5, loc= "upper right")
plt.ylim(-250, 1500)
plt.xlabel("time steps",fontsize=14)
plt.ylabel("proteins per cell ± s.e.m.",fontsize=14)
plt.title("average of the population",fontsize=14)
plt.show()

#PLOT+- SD
plt.plot(time_to_plot,proteins[0],label="lacI", color="tomato",lw=3)
plt.fill_between(time_to_plot,np.asarray(proteins[0])-np.asarray(sd[0]) , np.asarray(proteins[0])+ np.asarray(sd[0]), color="tomato", alpha=0.2, edgecolor='none')
plt.plot(time_to_plot,proteins[1],label="tetR", color="c",lw=3)
plt.fill_between(time_to_plot,np.asarray(proteins[1])-np.asarray(sd[1]) , np.asarray(proteins[1])+ np.asarray(sd[1]), color="c", alpha=0.2, edgecolor='none')
plt.plot(time_to_plot,proteins[2],label="cI", color="y",lw=3)
plt.fill_between(time_to_plot,np.asarray(proteins[2])-np.asarray(sd[2]) , np.asarray(proteins[2])+ np.asarray(sd[2]), color="y", alpha=0.2, edgecolor='none')
plt.legend(fontsize=12,framealpha=0.5, loc= "upper right")
plt.ylim(-250, 1500)
plt.xlabel("time steps",fontsize=14)
plt.ylabel("proteins per cell ± s.d.",fontsize=14)
plt.title("Average of the population",fontsize=14)
plt.show()


for i,name in enumerate(dict_to_plot):
    if i%3 == 0:
        plt.plot(time_to_plot, dict_to_plot[name], color="tomato")
    elif i%3 == 1:
        plt.plot(time_to_plot, dict_to_plot[name], color="c")
    elif i%3 == 2:
        plt.plot(time_to_plot, dict_to_plot[name], color="y")
plt.xlabel("time steps", fontsize=14)
plt.ylabel("proteins per cell", fontsize=14)
plt.title("Signal of 10 individual cells", fontsize=14)
plt.legend(fontsize=12,framealpha=0.5)
plt.ylim(-250, 6200)
plt.show()



#################################
#  PLOTING THE AUTOCORRELATION  #
#################################


#function from: https://stackoverflow.com/questions/14297012/estimate-autocorrelation-using-python
def acorr(x):
    n = len(x)
    var = variance(x)
    x = x-mean(x)
    r = np.correlate(x, x, mode = 'full')[-n:] #we take only one part of the results
    result = r/(var*(np.arange(n, 0, -1))) 
    return result


store_acorr_lacI=[]
store_acorr_tetR=[]
store_acorr_cI=[]
for i in dict_to_plot:
    if str(i)[0:2]=="la": #here we decide which of the three proteins we want to plot
        store_acorr_lacI.append(acorr(dict_to_plot[i]))
    if str(i)[0:2]=="te": #here we decide which of the three proteins we want to plot
        store_acorr_tetR.append(acorr(dict_to_plot[i]))
    if str(i)[0:2]=="cI": #here we decide which of the three proteins we want to plot
        store_acorr_cI.append(acorr(dict_to_plot[i]))

average_acorr_lacI=[np.mean(k) for k in zip(*store_acorr_lacI)]
average_acorr_tetR=[np.mean(k) for k in zip(*store_acorr_tetR)]
average_acorr_cI=[np.mean(k) for k in zip(*store_acorr_cI)]
# To not have messy results at the end I generate data for more time steps, like 200, but I only plot the first 100
plt.ylim(-0.5, 1.1)
plt.plot(range(100),average_acorr_lacI[:100], color="tomato", label="lacI",lw=3)
plt.plot(range(100),average_acorr_tetR[:100], color="c", label="tetR",lw=3)
plt.plot(range(100),average_acorr_cI[:100], color="y", label="cI",lw=3)
plt.xlabel( "shifted time steps",fontsize=14)
plt.ylabel('autocorrelation',fontsize=14)
plt.title("Averaged autocorrelation (200 cells)", fontsize=14)
plt.legend(fontsize=12, loc="upper right")
plt.show()

#######################
#  TRACKING DIVISION  #
#######################

division_track=[]
for cell in range(len(grid)):
    division_track.append(grid[cell]["div"]/time_limit)
plt.hist(division_track)
#values, counts = np.unique(division_track, return_counts=True)
#counts_freq=[x/len(division_track) for x in counts]

#plt.vlines(values, 0, counts_freq, color='C0', lw=20)
plt.ylabel("frequency")
# optionally set y-axis up nicely
#plt.ylim(0, max(counts_freq) * 1.06)
plt.show()

#Now we check if that is normal distribution
shapiro(division_track)






