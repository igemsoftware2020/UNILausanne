# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 09:05:55 2020

@author: pablo
"""

#NOTE FROM THE FUTURE, DO NOT ASSING THE COLOR BASED ON THE HIGHEST PROTEIN
#changes from the original version (200812)
#I simplify cell division just to 1 step (intead of first activation and later division)
#Cell division: the cell will chose a random position around it (let´s see how many layers around), and if it is free it will divide, otherwise it wont. As in "How spatial structure and neighbor uncertainty promote mutualists and weaken black queen effects"
#now cells cannot move without division
#Simplify the repressilator by considering only the proteins
#could finding of this model, the more layers away you can go when dividing, the thiker rings you get 


#I import the packages I will use
import numpy as np
import random
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import colors

#################
# Repressilator #
#################


def repressilator(z, t):
    p_lacI=z[0]
    p_tetR=z[1]
    p_cI=z[2]
    #
    alpha,n, nIPTG, Kd=(20,2.4,2, 10^-10)
    dp_lacIdt = (alpha/(1+p_cI**n))*(1 - IPTG_0**nIPTG/(Kd**nIPTG + IPTG_0**nIPTG))
    
    #
    dp_tetRdt = alpha/(1+p_lacI**n)
    #
    dp_cIdt=alpha/(1+p_tetR**n)
    return[dp_lacIdt,dp_tetRdt,dp_cIdt]


##########################
# Individual Based Model #
##########################

## PARAMETERS ##
layers=2#(layers of positions where the cell can move when dividin, 1 layer would mean the 8 adjacent positions, 2 layers the 24 adjacent positions and so on)
IPTG_0=10
death_rate=0.3
grid_size=100 #here I set the grid size, which will be the number of rows and columsn of a square grid
time_steps=2001; time_step=0#numer of rounds and the initial round for the while loop
cell_0={"div":1, "rep":[0,10,0]} #in the repressilator we include [mlacI,placI,mtetR,ptetR,mcI,pcI], if we have synchronised them with IPTG we will 



## MATRICES WE NEED ##

#The grid where cells will lie:
grid = [ [ 0 for i in range(grid_size) ] for j in range(grid_size) ]

#the idea is to implement a nutrient gradient coming from above in the future, but for now:
#nutrient= np.full((grid_size, grid_size), 0.4) #a matrix with the amount of nutrients in each of the positions
#for now they are not consuming nutrients, improve this


## FUNCTIONS TO USE LATER ON ###

#First we need to seed the initial cells. We randomly place the inputed bacteria in the bottom row of the matrix
#---------------------------------------------------------------------------------------------------------------
def initiation (number):
    var=0
    while var< number:
        
        column= random.randint(0,grid_size-1)
        if grid[grid_size-1][column]==0: #CHANGE 50 BY COLUMN
            grid[grid_size-1][column]=cell_0.copy() #CHANGE 50 BY COLUMN
            var+=1
        else:
            pass
        

#Now we define how cells will DIVIDE:
#------------------------------------
        #division will be something like this: (Based on https://doi.org/10.1016/j.bej.2019.107305)
        # CELL + NUTRIENTS -> ACTIVATED CELL -> 2 CELLS

#So then I will need 3 functions for division:
        #activation of cells
        #division of cells
        #movement of cells, so they can shove the cells surrounding them when dividing.


    


#I create a list of directions that I will use for the division
direction=[[-1,0],[0,1],[1,0],[0,-1]]

#I create a function to move cells around one dividing cell.When dividing the new cell will first be able to occupy one of the immediate 8 adjacent spots if it is free, adn otherwise move ONLY 1 layer of cells around to make some space. This is just a function to implement later on in the division function.
#In the model we will consider that cells cannot pass from bottom to top (which would mean to go from the gut wall to the gut lumen), but they can move from left ot right(they are not constrained horizontally)
#So that, rows can only be in the range of the matrix, but for columns they just divide them by %grid_size, this would mean taking the remainder of the division, so if it is column 20 it will go back to column 0, 21 to column 1 and so on...






def division(row, column):
    if random.uniform(0,1)<grid[row][column]["div"]: #only excited cells (type 2)
        a=random.choice(range(-layers, layers+1)) #layers +1 because range doesn´t include the upper number
        b=random.choice(range(-layers, layers+1))
        if (row+a) in range(grid_size):
            if grid[row+a][(column+b)%grid_size]==0:
                    grid[row+a][(column+b)%grid_size]=grid[row][column].copy()
                    new_rep=[x/2 for x in grid[row][column]["rep"]]
                    variability=np.random.normal(1, 0.2, 1)[0]
                    new_rep_1=[x*variability for x in new_rep]
                    new_rep_2=[x*(2-variability) for x in new_rep]
                    grid[row+a][(column+b)%grid_size]["rep"]=new_rep_2
                    grid[row][column]["rep"]=new_rep_1


#another function for cell death (just set the cell and all its parameters to 0)
#-------------------------------
def death(row, column):
    grid[row][column]=0

    
    

#How to plot the matrix
#----------------------
    #Even though we have cells as dictionaries, we need to just have some values to be able to plot the grid properly, what we do for now is basically check which is the higher protien of the three in the repressilator of a particular cell, and based on that assign 1, 2 or 3 to that particular position in the grid (only for plotting)
def plotable(grid_0):
    grid_1 = [row[:] for row in grid_0]
    for r in range(grid_size):
        for c in range(grid_size):
            if grid_1[r][c] != 0:
                prot=[grid_1[r][c]["rep"][0],grid_1[r][c]["rep"][1],grid_1[r][c]["rep"][2]]
                grid_1[r][c]=1+prot.index(max(prot))
    return grid_1
    



#################################
## NOW WE START THE SIMULATION ##
#################################

#we initialte the grid with the number of cells we want:
initiation(1)


#we plot the grid before starting (just to see where the cells went randomly)
plt.imshow(plotable(grid),cmap = colors.ListedColormap(["black",'tomato','turquoise',"y"]))
plt.title('before starting')
plt.colorbar()
plt.show()


time_step=0 #we also set this above, I needed it here hen removing IPTG

#Now perform the simulation, basically at each time step, each bacteria can die, run a bit the repressialtor, be activated or divide
while time_step<time_steps: #the simulation will run until the number of rounds previously set
    grid_copy=np.copy(grid)#I copy the matrix to not loop and modify at the same time (this should be done in other way I think)
    for i in range(grid_size):
        for z in range(grid_size):
            if grid_copy[i][z]!= 0:
                if random.uniform(0, 1)<death_rate: #It will die if a random number is lower than the death rate we set (+ toxin presence in that position)
                    death(i,z)
                else:
                    wi = odeint(repressilator, grid[i][z]["rep"], [0,0.02] )
                    grid[i][z]["rep"]=wi[1]
                    division(i,z)


#now we just plot the grid at each time stpel.
    
    if time_step%5==0: #I print the matrix fro every 100, 200, 300....round
        plt.imshow(plotable(grid), cmap = colors.ListedColormap(["black",'tomato','turquoise',"y"]),vmin=0,vmax=3) #lacI==RED; tetR==BLUE; cI==YELLOW
        plt.title("time: " + format(time_step))
        plt.colorbar()
        #plt.savefig("plots/fig__" + format(time_step)+".png")
        plt.show()
    
    else:
        pass
    time_step+=1

#In case we want to remove IPTG from each cell: 
#for x in range(grid_size):
    #for y in range(grid_size):
        #if grid[x][y]!=0:
            #grid[x][y]["rep"][6]=0