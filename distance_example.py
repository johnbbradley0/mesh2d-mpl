
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
from scipy import linspace, array, meshgrid, where
from matplotlib.pyplot import figure, subplot, contourf, colorbar, close, title

# my imports 
import distance 
BoxDistance, DiskDistance = distance.BoxDistance, distance.DiskDistance


# Global variables
# ================================================================================
dtotal = BoxDistance([-2, 2, 0, 2])
dcloud = BoxDistance([-.5, .5, .2, .4]) + BoxDistance([-1.5, -.45, .2, 1.4]) 
xs = ([.5, .7,  .9, 1.1, 1.3, 1.5])
rs = [.18, .22, .25, .25, .22, .18]
for x, r in zip(xs, rs):
    dcloud += DiskDistance([x,1.25],r)

dclear = dtotal - dcloud



# Example
# ================================================================================
close(0)
x_ = linspace(-2.5, 2.5, 400)
y_ = linspace(-0.5, 2.5, 400)
x, y = meshgrid(x_,y_)
f0 = figure(0, (8,14), facecolor='w')

s1 = subplot(311); 
dtotal.plot(x,y,"Total")

s2 = subplot(312)
dcloud.plot(x,y,"Clouds")

s3 = subplot(313); 
dclear.plot(x,y,"Clear")

f0.show()



