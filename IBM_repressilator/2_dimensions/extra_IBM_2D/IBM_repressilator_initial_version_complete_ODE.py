# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:39:37 2020

@author: pablo
"""

#One of the initial versions of our IBM in 2 dimension
#The ODE we used to update the proteins every time step were more complete (ELowitz et Leiber 2000), but also because of that was more costly
#It has some probles that we later realised such as the fact that cells are coloured by the higher protein and not by the intensity that the amount of protein will give them


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
   
    mRNA_lacI = z[0]
    p_lacI    = z[1]
    mRNA_tetR = z[2]
    p_tetR    = z[3]
    mRNA_cI   = z[4]
    p_cI      = z[5]
    IPTG      = z[6]
    
    alpha0 = 0.216
    n      = 2
    nIPTG  = 1
    beta   = 1
    Kd     = 10^-10
    
    #IPTG
    dIPTG        = IPTG * (1 - p_lacI**nIPTG/(Kd**nIPTG + p_lacI**nIPTG) )
    
    #LacI
    alpha_lacI   = 100
    dmRNA_lacIdt = alpha_lacI/(1 + p_cI**n) + alpha0 - mRNA_lacI
    dp_lacIdt    = - beta*( p_lacI - mRNA_lacI*(1 - IPTG**nIPTG/(Kd**nIPTG + IPTG**nIPTG)))
    
    #TetR
    alpha_tetR   = 100
    dmRNA_tetRdt = alpha_tetR/(1 + p_lacI**n) + alpha0 - mRNA_tetR
    dp_tetRdt    = - beta*(p_tetR -mRNA_tetR)

    #cI
    alpha_cI     = 100
    dmRNA_cIdt   = alpha_cI/(1+(p_tetR)**n)+ alpha0 - mRNA_cI
    dp_cIdt      = - beta*(p_cI-mRNA_cI)
    
    
    return[dmRNA_lacIdt, dp_lacIdt, dmRNA_tetRdt, dp_tetRdt, dmRNA_cIdt, 
           dp_cIdt, dIPTG]



##########################
# Individual Based Model #
##########################

## PARAMETERS ##
IPTG_0=0
death_rate=0.005
grid_size=30 #here I set the grid size, which will be the number of rows and columsn of a square grid
time_steps=101; time_step=0#numer of rounds and the initial round for the while loop
nutrients_newborn=0.05#nutrients that consumes a newborncell 
cell_0={"type":1, "react":0.5, "div":0.5, "rep":[0,0,10,10,0,0,IPTG_0]} #in the repressilator we include [mlacI,placI,mtetR,ptetR,mcI,pcI], if we have synchronised them with IPTG we will 



## MATRICES WE NEED ##

#The grid where cells will lie:
grid = [ [ 0 for i in range(grid_size) ] for j in range(grid_size) ]

#the idea is to implement a nutrient gradient coming from above in the future, but for now:
nutrient= np.full((grid_size, grid_size), 0.4) #a matrix with the amount of nutrients in each of the positions



## FUNCTIONS TO USE LATER ON ###

#First we need to seed the initial cells. We randomly place the inputed bacteria in the bottom row of the matrix
#---------------------------------------------------------------------------------------------------------------
def initiation (number):
    var=0
    while var< number:
        
        column= random.randint(0,grid_size-1)
        if grid[grid_size-1][column]==0:
            grid[grid_size-1][column]=cell_0.copy()
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


def activation( r, c):
    if grid[r][c]["type"]==1:
        if random.uniform(0,1)<grid[r][c]["react"]*nutrient[r][c]:
            if nutrient[r][c]>=nutrients_newborn:
                grid[r][c]["type"]+=1  
                nutrient[r][c]-=nutrients_newborn
    else:
        pass


#I create a list of directions that I will use for the division
direction=[[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]]

#I create a function to move cells around one dividing cell.When dividing the new cell will first be able to occupy one of the immediate 8 adjacent spots if it is free, adn otherwise move ONLY 1 layer of cells around to make some space. This is just a function to implement later on in the division function.
#In the model we will consider that cells cannot pass from bottom to top (which would mean to go from the gut wall to the gut lumen), but they can move from left ot right(they are not constrained horizontally)
#So that, rows can only be in the range of the matrix, but for columns they just divide them by %grid_size, this would mean taking the remainder of the division, so if it is column 20 it will go back to column 0, 21 to column 1 and so on...


def move(row,column):
    for dir_ in direction:
        row_2=(row+dir_[0])
        column_2=(column+dir_[1])
        if row_2 in range(grid_size) and column_2 in range (grid_size):
            if row_2+dir_[0] in range(grid_size):
                if grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]==0:
                    grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]=grid[row_2][column_2]
                
                break



def division(row, column):
    if grid[row][column]["type"]==2 and random.uniform(0,1)<grid[row][column]["div"]: #only excited cells (type 2)
        random.shuffle(direction)
        for dir_ in direction:
            if (row+dir_[0]) in range(grid_size) and (column+dir_[1]) in range(grid_size):
                if grid[(row+dir_[0])][(column+dir_[1])]==0:
                    grid[row][column]["type"]=1
                    grid[(row+dir_[0])][(column+dir_[1])]=grid[row][column].copy()
                    new_rep=[x/2 for x in grid[row][column]["rep"]]
                    variability=np.random.normal(1, 0.2, 1)[0]
                    new_rep_1=[x*variability for x in new_rep]
                    new_rep_1[6]=IPTG_0#somehow the new cell will have the same ammount of IPTG as outside
                    new_rep_2=[x*(2-variability) for x in new_rep]
                    new_rep_2[6]=IPTG_0#somehow the new cell will have the outside concentration of IPTG
                    grid[(row+dir_[0])][(column+dir_[1])]["rep"]=new_rep_2
                    grid[row][column]["rep"]=new_rep_1
                    
                    break
                else:
                    move(row,column)
                    break

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
                prot=[grid_1[r][c]["rep"][1],grid_1[r][c]["rep"][3],grid_1[r][c]["rep"][5]]
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

#I am not sure whether time for integration sould be like that, but for now this is working
t=np.linspace(0,time_steps/10,time_steps+1)
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
                    wi = odeint(repressilator, grid[i][z]["rep"], [t[time_step],t[time_step+1]] )
                    grid[i][z]["rep"]=wi[1]
                    activation(i,z)
                    division(i,z)

#now we just plot the grid at each time stpel.
    
    if time_step%1==0: #I print the matrix fro every 100, 200, 300....round
        plt.imshow(plotable(grid), cmap = colors.ListedColormap(["black",'tomato','turquoise',"y"]),vmin=0,vmax=3) #lacI==RED; tetR==BLUE; cI==YELLOW
        plt.title("time: " + format(time_step))
        plt.colorbar()
        plt.show()
        
    
    else:
        pass
    time_step+=1

#In case we want to remove IPTG from each cell: 
#for x in range(grid_size):
    #for y in range(grid_size):
        #if grid[x][y]!=0:
            #grid[x][y]["rep"][6]=0