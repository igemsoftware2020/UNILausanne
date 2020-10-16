# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 09:05:55 2020

@author: pablo
"""

#here I dont run the repressilator at each step to simplify, in the one from yesterday I have everything

#I import the packages I will use
import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import shapiro
from matplotlib import colors



##########################
# Individual Based Model #
##########################

normality_check=[]

## PARAMETERS ##
grid_size=100 #here I set the grid size, which will be the number of rows and columsn of a square grid
time_steps=501; time_step=0#numer of rounds and the initial round for the while loop
cell_0={"div":0}
                 

modes=[0.4,0.8,0.24,1,2.01,2.05,2.075,2.1,3.0101,3.0301,3.0501,3.0701,3.1001,3.0105,3.0305,3.0505,3.0705,3.1005]
death_rates=[round(x*0.05,2) for x in range(20)]
div_probs=[0.85,0.5,0.15]
initial_cells=[1,10,100]

for mode in modes:
    for death_rate in death_rates:
        for mean_prob in div_probs:
            for in_cells in initial_cells:
                
                
                
                
                ## MATRICES WE NEED ##
                
                #The grid where cells will lie:
                grid = [ [ 0 for i in range(grid_size) ] for j in range(grid_size) ]
                
                #the idea is to implement a nutrient gradient coming from above in the future, but for now:
                #nutrient= np.full((grid_size, grid_size), 0.4) #a matrix with the amount of nutrients in each of the positions
                #for now they are not consuming nutrients, I think they might not be important in this simulation
                
                
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
                        
                
                
                #I create a list of directions that I will use for the division MODE 0
                direction_4= [[-1,0],[0,1],[1,0],[0,-1]]
                direction_8=[[-1,-1],[-1,0],[-1,1],[0,1],[1,1],[1,0],[1,-1],[0,-1]]
                direction_24=[[-2,-2],[-1,-2],[0,-2],[1,-2],[2,-2],[-2,-1],[-1,-1],[0,-1],[1,-1],[2,-1],[-2,0],[-1,0],[1,0],[2,0],[-2,1],[-1,1],[0,1],[1,1],[2,1],[-2,2],[-1,2],[0,2],[1,2],[2,2]]
                
                #I will use the shove function for division MODE 1
                def shove(row,column):
                    for dir_ in direction_8:
                        row_2=(row+dir_[0])
                        column_2=((column+dir_[1])%grid_size)
                        if row_2 in range(grid_size) and row_2+dir_[0] in range(grid_size):
                            if grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]==0 and grid[row_2][column_2]!=0:#theoreically I shouldn need to ask the second condition,but was not working I dont know why´
                                grid[(row_2+dir_[0])][(column_2+dir_[1])%grid_size]=grid[row_2][column_2].copy()
                                grid[row_2][column_2]=grid[row][column].copy()
                            
                                break
                
                #In the model we will consider that cells cannot pass from bottom to top (which would mean to go from the gut wall to the gut lumen), but they can move from left ot right(they are not constrained horizontally)
                #So that, rows can only be in the range of the matrix, but for columns they just divide them by %grid_size, this would mean taking the remainder of the division, so if it is column 20 it will go back to column 0, 21 to column 1 and so on...
                
                def move(row,column):
                    if random.uniform(0,1)<move_prob:
                        a=direction[random.choice(range(len(direction)))]
                        if (row+a[0]) in range(grid_size):
                            if grid[(row+a[0])][(column+a[1])%grid_size]==0:
                                grid[(row+a[0])][(column+a[1])%grid_size]=grid[row][column].copy()
                                grid[row][column] = 0
                #DESCRIPTION OF DIVISION OPTIONS
                
                #MODE=0, we will be looping in all those positions until we find an empty space
                    # Here we have to chose direction 4,8, or 24
                direction=direction_4
                if mode == 0.4:
                    direction= direction_4
                elif mode == 0.8:
                    direction= direction_8
                elif mode == 0.24:
                    direction= direction_24
                #MODE=1, can divide to the 8 adjacent and shove one extra layer (they can shove in the same direction)
                
                #MODE=2, cells can move to one of the four adjacent positions, but can move without division. Set the move prob here
                #here we need to provide the moving probability
                
                elif mode == 2.01:
                    move_prob = 0.1
                elif mode == 2.05:
                    move_prob = 0.5
                elif mode == 2.075:
                    move_prob = 0.75
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
                elif mode ==3.0701:
                    layers=7
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
                
                
                def division(mode,row, column):

                    div_prob=np.random.normal(mean_prob, 0.1, 1)[0] #I want to try what sara suggested, lets keep sd=0.1 for now

                    if mode==0:
                        if random.uniform(0,1)<div_prob: 
                            random.shuffle(direction)
                            for dir_ in direction:
                                if (row+dir_[0]) in range(grid_size):
                                    if grid[(row+dir_[0])%grid_size][(column+dir_[1])%grid_size]==0:
                                        grid[row][column]["div"] += 1
                                        grid[row+dir_[0]][(column+dir_[1])%grid_size]=grid[row][column].copy()
                                        
                                        break
                        
                      
                                   
                    elif mode == 1:
                        if random.uniform(0,1)< div_prob: 
                            random.shuffle(direction_8)
                            could_divide=0
                            for dir_ in direction_8:
                                if (row+dir_[0]) in range(grid_size) :
                                    if grid[(row+dir_[0])%grid_size][(column+dir_[1])%grid_size]==0:
                                        grid[row][column]["div"] += 1
                                        grid[row+dir_[0]][(column+dir_[1])%grid_size]=grid[row][column].copy()
                                        
                                        could_divide +=1
                                        break
                
                                    if could_divide == 0:
                                        shove(row,column)
                                        break
                    elif mode == 2:
                        if random.uniform(0,1)<div_prob: 
                            a=direction[random.choice(range(len(direction_4)))]
                            if (row+a[0]) in range(grid_size):
                                if grid[(row+a[0])][(column+a[1])%grid_size]==0:
                                    grid[row][column]["div"] += 1
                                    grid[(row+a[0])][(column+a[1])%grid_size]=grid[row][column].copy()
                                 
                        move(row,column)
                    
                    if mode == 3:
                        if random.uniform(0,1)<div_prob: 
                            attempt=0
                            while attempt<attempts:
                                a=random.choice(range(-layers, layers+1)) #layers +1 because range doesn´t include the upper number
                                b=random.choice(range(-layers, layers+1))
                                if (row+a) in range(grid_size):
                                    if grid[row+a][(column+b)%grid_size]==0:
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
                
                def plotable_cell(grid_0):
                    grid_1 = [row[:] for row in grid_0]
                    for r in range(grid_size):
                        for c in range(grid_size):
                            if grid_1[r][c] != 0:
                                grid_1[r][c]=1
                    return grid_1
                
                def plotable_div(grid_0):
                    grid_1 = [row[:] for row in grid_0] #we copy the matrix
                    for r in range(grid_size):
                        for c in range(grid_size):
                            if grid_1[r][c] != 0:
                                grid_1[r][c]=grid_1[r][c]["div"]
                    return grid_1
                    
                
                #################################
                ## NOW WE START THE SIMULATION ##
                #################################
                
                #we initialte the grid with the number of cells we want:
                initiation(in_cells)
                
                
                
                
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
                                    division(int(mode),i,z)
                
                
                #now we just plot the grid at each time stpel.
                    if time_step==50 or time_step%100==0: #I  print the matrix fro every 100, 200, 300....round
                        plt.imshow(plotable_div(grid))
                        plt.title("time: " + format(time_step))
                        plt.colorbar()
                        plt.savefig("201007_plots/fig_" + format(mode) + "_"+ format(mean_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_"+ format(time_step)+".png")
                        plt.show()
                 
                        
                        division_track=[]
                        for thing in grid:
                            for cell in thing:
                                if cell != 0:
                                    try: #still geting some 0 idkwhy so I added this try except to solve it for now
                                        division_track.append(cell["div"]/time_steps)
                                    except:
                                        continue
                        lets_see=0
                        for i in division_track:
                            if i:
                                lets_see +=1
                        if lets_see>=10:
                            
                            plt.hist(division_track)
                            plt.ylabel("frequency")
                            plt.savefig("201007_plots/hist_" + format(mode) + "_"+ format(mean_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_" + format(time_step)+".png")
                            plt.show
                            
                            pd.DataFrame(division_track).to_csv("201007_data/"+ format(mode) + "_"+ format(mean_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_" + format(time_step) +".csv")
                            #dt=np.array(division_track)
                            #np.savez("data_2/"+ format(mode) + "_"+ format(mean_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_" + format(time_step), dt)
                            #a = np.load("Myfile.npz") #make sure you use the .npz!
                            #b = a['arr_0']
                            
                           
                            test = shapiro(division_track)
                            normality_check.append("SW_test_"+ format(mode) + "_"+ format(mean_prob) + "_" +format(death_rate)+"_" + format(in_cells)+ "_" + format(time_step) + ",stat," + str(test[0]) + ",pvalue," + str(test[1]))
                            
                        
                    time_step+=1
                
                
                
                #####################################
                # LET´S TRACK HOW MUCH CELLS DIVIDE #
                #####################################
pd.DataFrame(normality_check).to_csv("201007_data/normality_check.csv")
        
