# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 09:05:55 2020

@author: pablo
"""

#Model to explore different options for cell division in our IBM
#We have decribed the 5 modes of division in our wiki, check that out in our wiki (link in the readme of our repo)


#I import the packages I will use (and some more)
import numpy as np
import random
import matplotlib.pyplot as plt
import scipy
from scipy.integrate import odeint
from matplotlib import colors
import pandas as pd
from scipy.stats import shapiro
import statsmodels.api as sm 
import pylab as py 

#################
# Repressilator #
#################


def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    #
    alpha,n, nIPTG, Kd=(216,2.4, 1, 1e-10) #Parameters from Elowitz et al 2000.Except for Kd that is taken fro Daber et al 2007
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG)) #In case we want to coordinate with IPTG
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n)
    #
    dp_cIdt=alpha/(1+p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


##########################
# Individual Based Model #
##########################



## PARAMETERS ##
#--------------------------------------------------------------------------------------------------------------

IPTG_0=0
death_rate=0
grid_size=100 #here I set the grid size, which will be the number of rows and columsn of a square grid
time_steps=91; time_step=0#numer of rounds and the initial round for the while loop

cell_0={"div":0, "rep":[0,10,0]} #in the repressilator we include [LacI,TetR,CI], we make up the initial values

## FUNCTIONS ##
#--------------------------------------------------------------------------------------------------------------

#Function to initialise the matrix trandomly in the row we want with the desited number of bacteria, now I changed it to just put one bacterium in the center
def initiation (number):
    var=0
    while var< number:
        
        column= random.randint(0,grid_size-1)
        if grid[grid_size/2][grid_size/2]==0: #CHANGE griz_size/2 BY COLUMN AND SELECT THE ROW
            grid[grid_size/2][grid_size/2]=cell_0.copy() #CHANGE grid_size/2 BY COLUMN AND SELECT THE ROW
            var+=1
        else:
            pass
        


#I create a list of directions that I will use for the division MODE 0
direction_4= [[-1,0],[0,1],[1,0],[0,-1]]
direction_8=[[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]]
direction_24=[[-2,-2],[-1,-2],[0,-2],[1,-2],[2,-2],[-2,-1],[-1,-1],[0,-1],[1,-1],[2,-1],[-2,0],[-1,0],[1,0],[2,0],[-2,1],[-1,1],[0,1],[1,1],[2,1],[-2,2],[-1,2],[0,2],[1,2],[2,2]]

#I will use the shove function for division MODE 1
#I create a function to move cells around one dividing cell.When dividing the new cell will first be able to occupy one of the immediate 8 adjacent spots if it is free, adn otherwise move ONLY 1 layer of cells around to make some space. This is just a function to implement later on in the division function.
def shove(row,column):
    for dir_ in direction_8:
        row_2=(row+dir_[0])
        column_2=((column+dir_[1])%grid_size)
        if row_2 in range(grid_size) and row_2+dir_[0] in range(grid_size):
            if grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]==0 and grid[row_2][column_2]!=0:#theoreically I shouldn need to ask the second condition,but was not working I dont know why´
                grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]=grid[row_2][column_2].copy()
                grid[row_2][column_2]=grid[row][column].copy()
                new_rep=[x/2 for x in grid[row][column]["rep"]]
                variability=np.random.normal(1, 0.2, 1)[0]
                new_rep_1=[x*variability for x in new_rep]
                new_rep_2=[x*(2-variability) for x in new_rep]
                grid[row_2][column_2]["rep"]=new_rep_2
                grid[row][column]["rep"]=new_rep_1
                break

def move(row,column):
    if random.uniform(0,1)<move_prob:
        a=direction_4[random.choice(range(len(direction_4)))]
        if (row+a[0]) in range(grid_size):
            if grid[(row+a[0])][(column+a[1])%grid_size]==0:
                grid[(row+a[0])][(column+a[1])%grid_size]=grid[row][column].copy()
                grid[row][column] = 0


#In the model we will consider that cells cannot pass from bottom to top (which would mean to go from the gut wall to the gut lumen), but they can move from left ot right(they are not constrained horizontally)
#So that, rows can only be in the range of the matrix, but for columns they just divide them by %grid_size, this would mean taking the remainder of the division, so if it is column 20 it will go back to column 0, 21 to column 1 and so on...
#Here we have plenty modes of division
def division(mode,row, column):
    variability=0.5# np.random.beta(50,50)
    if mode==0:
        if random.uniform(0,1)<div_prob: #only excited cells (type 2)
            random.shuffle(direction)
            for dir_ in direction:
                if (row+dir_[0]) in range(grid_size):
                    if grid[(row+dir_[0])%grid_size][(column+dir_[1])%grid_size]==0:
                        grid[row][column]["div"] += 1
                        grid[row+dir_[0]][(column+dir_[1])%grid_size]=grid[row][column].copy()
                        new_rep=grid[row][column]["rep"].copy()
                        new_rep_1=[x*variability for x in new_rep]
                        new_rep_2=[x*(1-variability) for x in new_rep]
                        grid[row+dir_[0]][(column+dir_[1])%grid_size]["rep"]=new_rep_2
                        grid[row][column]["rep"]=new_rep_1
                        break
    elif mode == 1:
        if random.uniform(0,1)< div_prob: #only excited cells (type 2)
            random.shuffle(direction_8)
            could_divide=0
            for dir_ in direction_8:
                if (row+dir_[0]) in range(grid_size) :
                    if grid[(row+dir_[0])%grid_size][(column+dir_[1])%grid_size]==0:
                        grid[row][column]["div"] += 1
                        grid[row+dir_[0]][(column+dir_[1])%grid_size]=grid[row][column].copy()
                        new_rep=grid[row][column]["rep"].copy()
                        new_rep_1=[x*variability for x in new_rep]
                        new_rep_2=[x*(1-variability) for x in new_rep]
                        grid[row+dir_[0]][(column+dir_[1])%grid_size]["rep"]=new_rep_2
                        grid[row][column]["rep"]=new_rep_1
                        could_divide +=1
                        break
    
                    if could_divide == 0:
                        shove(row,column)
                        break
    elif mode == 2:
        if random.uniform(0,1)<div_prob: #only excited cells (type 2)
            a=direction_4[random.choice(range(len(direction_4)))]
            if (row+a[0]) in range(grid_size):
                if grid[(row+a[0])][(column+a[1])%grid_size]==0:
                    grid[row][column]["div"] += 1
                    grid[(row+a[0])][(column+a[1])%grid_size]=grid[row][column].copy()
                    new_rep=grid[row][column]["rep"].copy()
                    new_rep_1=[x*variability for x in new_rep]
                    new_rep_2=[x*(1-variability) for x in new_rep]
                    grid[(row+a[0])][(column+a[1])%grid_size]["rep"]=new_rep_2
                    grid[row][column]["rep"]=new_rep_1.copy()
            
                
                (row,column)
    
    
    if mode == 3:
            if random.uniform(0,1)<div_prob: 
                attempt=0
                while attempt<attempts:
                    a=random.choice(range(-layers, layers+1)) #layers +1 because range doesn´t include the upper number
                    b=random.choice(range(-layers, layers+1))
                    if (row+a) in range(grid_size):
                        if grid[row+a][(column+b)%grid_size]==0:
                            grid[row][column]["div"] += 1
                            grid[row+a][(column+b)%grid_size]=grid[row][column].copy()
                            new_rep=grid[row][column]["rep"].copy()
                            new_rep_1=[x*variability for x in new_rep]
                            new_rep_2=[x*(1-variability) for x in new_rep]
                            grid[row+a][(column+b)%grid_size]["rep"]=new_rep_2
                            grid[row][column]["rep"]=new_rep_1
                            grid[row+a][(column+b)%grid_size]=grid[row][column].copy()
                          
                            break
                        else:
                            attempt +=1
                       
            
            
#another function for cell death (just set the cell and all its parameters to 0)
#-------------------------------
def death(row, column):
    grid[row][column]=0

    
    

#How to plot the matrix
#----------------------
    #Even though we have cells as dictionaries, we need to just have some values to be able to plot the grid properly, what we do for now is basically check which is the higher protien of the three in the repressilator of a particular cell, and based on that assign 1, 2 or 3 to that particular position in the grid (only for plotting)
def plotable_lacI(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                grid_1[r][c]=grid_1[r][c]["rep"][0]*40
    return grid_1
    
def plotable_tetR(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                grid_1[r][c]=grid_1[r][c]["rep"][1]*40
    return grid_1

def plotable_cI(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                grid_1[r][c]=grid_1[r][c]["rep"][2]*40
    return grid_1

def plotable_div(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                grid_1[r][c]=grid_1[r][c]["div"]
    return grid_1



def plotable(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                grid_1[r][c]=1
    return grid_1
#--------------------------------------------------------------------------------------------------------------


#All the options we want to try
modes=[0.4,0.8,0.24,1,2.05,2.1,3.0101,3.0301,3.0501,3.1001]
in_cells=1
move_=[0]
death_rates=[0,0.1,0.2,0.3]
div_probs=[1,0.8,0.5]


for mode in modes:
    for move_prob in move_:
        for death_rate in death_rates:
            for div_prob in div_probs:
                #We create the grid where cells will lie:
                grid = [ [ 0 for i in range(grid_size) ] for j in range(grid_size) ]
                
                #DESCRIPTION OF DIVISION OPTIONS
                #MODE=0, we will be looping in all those positions until we find an empty space
                        # Here we have to chose direction 4,8, or 24
        
                if mode == 0.4:
                    direction= direction_4
                elif mode == 0.8:
                    direction= direction_8
                elif mode == 0.24:
                    direction= direction_24
                #MODE=1, can divide to the 8 adjacent and shove one extra layer (they can shove in the same direction)
                
                #MODE=2, cells can move to one of the four adjacent positions, but can move without division. Set the move prob here
                #here we need to provide the moving probability

                elif mode == 2.05:
                    move_prob = 0.5

                elif mode == 2.1:
                    move_prob = 1
                
                
                #MODE=3, we select the number of layers where a cell can go when dividing and then it choses a position and if it is free it divides
                    #Here we need to provide the number of layers
                    #we can also have more than one attempt
                
                elif mode ==3.0101:
                    layers=1
                    attempts=1
                elif mode ==3.0301:
                    layers=3
                    attempts=1
                elif mode ==3.0501:
                    layers=5
                    attempts=1

                elif mode ==3.1001:
                    layers=10
                    attempts=1

                elif mode ==3.0105:
                    layers=1
                    attempts=5
                elif mode ==3.0305:
                    layers=3
                    attempts=5
                elif mode ==3.0505:
                    layers=5
                    attempts=5
                elif mode ==3.0705:
                    layers=7
                    attempts=5
                elif mode ==3.1005:
                    layers=10
                    attempts=5
                
                elif mode ==3.2005:
                    layers=20
                    attempts=5
                
                #################################
                ## NOW WE START THE SIMULATION ##
                #################################
                
                #we initialte the grid with the number of cells we want:
                initiation(in_cells) #we initialise the matrix
                time_step=0 #we set the time step to 0
                #Now perform the simulation, basically at each time step, each bacteria can die, run a bit the repressialtor, be activated or divide
                while time_step<time_steps: #the simulation will run until the number of rounds previously set
                    grid_copy=np.copy(grid)#I copy the matrix to not loop and modify at the same time (this should be done in other way I think)
                    #I create this to loop in random order! Important to not have undesired patterns.
                    size_1=list(range(grid_size))
                    size_2=list(range(grid_size))
                    random.shuffle(size_1)
                    random.shuffle(size_2)
                    for i in size_1:
                        for z in size_2:
                            if grid_copy[i][z]!= 0:
                                if random.uniform(0, 1)<death_rate: #It will die if a random number is lower than the death rate we set (+ toxin presence in that position)
                                    death(i,z)
                                else:
                                    #here we run repressialtor
                                    wi = odeint(repressilator, grid[i][z]["rep"], [0,0.7] )
                                    grid[i][z]["rep"]=wi[1]
                                    division(int(mode),i,z)
                                    move(i,z)
                    
                    if time_step%45==0: #I  print the matrix fro every 45, 90,....round, whenever the reminder equals 0!
                        #Here we plot the expression of one of the proteins, we select CI as has the same promoter as azurin
                        plt.imshow(plotable_cI(grid), cmap="magma")
                        plt.title("time: " + format(time_step))
                        plt.colorbar()
                        plt.savefig("201018/div" + format(mode) + "_"+ format(div_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_"+ format(time_step)+ "_"+format(move_prob)+".png", transparent=True)
                        plt.show()
                        
                        #Here we also plot where we have cells to see the shape of the colony
                        plt.imshow(plotable(grid), cmap="magma")
                        plt.title("time: " + format(time_step))
                        plt.colorbar()
                        plt.savefig("201018/colony" + format(mode) + "_"+ format(div_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_"+ format(time_step)+ "_"+format(move_prob)+".png", transparent=True)
                        plt.show()
                        
                    time_step+=1

