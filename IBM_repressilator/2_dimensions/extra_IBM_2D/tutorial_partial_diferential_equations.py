# -*- coding: utf-8 -*-
"""
Created on Sun Jul  5 19:34:12 2020

@author: pablo
"""

#At some point we discussed whether we should have diffusion of nutrients in our IBM 2D, for that we would use partial differential equations, so we tried this tutorila to have different amount of nutrients in each grid of the matrix. Bacteria would consume them, and meanwhile nutrients could difuse

#adapted from https://scipython.com/book/chapter-7-matplotlib/examples/the-two-dimensional-diffusion-equation/

import numpy as np
import matplotlib.pyplot as plt

# plate size, positions
w = h = 100

# intervals in x-, y- directions, positions
dx = dy = 1 #size of a grid position
# diffusivity of nutrients, I dont get which 
D = 10

#we set the minimal and maximal amount of nutrients
nut_min, nut_max = 0.4, 1

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = dx2 * dy2 / (2 * D * (dx2 + dy2))

#we initiatie the matrix
u0 = nut_min * np.ones((nx, ny))
u = u0.copy()

# Initial conditions - ring of inner radius r, width dr centred at (cx,cy) (mm)
# nutrients are coming from above, so, we will create a band of desired height
band_h=1 #positions

for i in range(int(band_h/dy)): # I add them in the top row
    for j in range(ny):
        u0[i,j] = nut_max
for i in range(int(h/2)): #I can also add them in half of the sides
    u0[i,0]=nut_max
    u0[i,-1]=nut_max

def do_timestep(u0, u):
    #If we want to regulate how much it difuses we can multiply it all by a constant somethign like 1.001
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] =  u0[1:-1, 1:-1] +  D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )
    #now I add this part to continusly rebuild the nutrients. 
    u[:int(band_h/dy),:]=nut_max
    for i in range(int(h/2)):
        u[i,0]=nut_max
        u[i,-1]=nut_max
    u0 = u.copy()
    return u0, u

# Number of timesteps
nsteps = 10001
# Output 4 figures at these timesteps
mfig = [0,1000,10000]
fignum = 0
fig = plt.figure()
for m in range(nsteps):
    u0, u = do_timestep(u0, u)
 
    if m in mfig:
        print(u)
        
        fignum += 1
        print(m, fignum)
        ax = fig.add_subplot(220 + fignum)
        im = ax.imshow(u.copy(), cmap=plt.get_cmap('hot'), vmin=nut_min,vmax=nut_max)
        ax.set_axis_off()
        ax.set_title('{:.1f} ms'.format(m*dt*1000))
fig.subplots_adjust(right=0.85)
cbar_ax = fig.add_axes([0.9, 0.15, 0.03, 0.7])
cbar_ax.set_xlabel('$T$ / K', labelpad=20)
fig.colorbar(im, cax=cbar_ax)
plt.show()
