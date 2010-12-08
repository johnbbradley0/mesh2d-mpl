
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
import scipy as sp
from enthought.mayavi import mlab


# import plane_parallel
from plane_parallel import Domain

# Global variables
# ================================================================================


# make a spatial_set object
a = 0.0  
b = 10.0
nz = 100

# function and class definitions 
# ================================================================================


        


# ====
# Example script
# ================================================================================

# imports
bgcolor = (0,0,0)
fgcolor = (1,1,1)
figsize = (1200, 900)


# make a spatial_set object
a = 0.0  
b = 10.0
nz = 100
hz = .6
hv = .35

# function for testing 
fun = lambda z: sp.exp(-(z-a))

domain = Domain(a,b,hz, hv)

mlab.close(0)    
f = mlab.figure(0, 
                bgcolor=bgcolor, 
                fgcolor=fgcolor, size=figsize)

# plot the figure
# mlab.clf(0)    

# plot a funciton in space
domain.plot_function(fun, tube_radius=.1)


# plot quivers
names = ['incoming', 'outgoing', 'domain']
colors = [(0,0,1), (1,0,0), (0,1,0)]
scale_factors = [.3, .3, .2]
for name, color, scale_factor in zip(names, colors, scale_factors):
    name += "_set_3d"
    domain.quiver3d(name, color=color, scale_factor=scale_factor)
    




















