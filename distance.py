
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
from scipy import linspace, array, meshgrid, where, sqrt
import numexpr as ne
import matplotlib.pyplot as plt

# my imports 


# Global variables
# ================================================================================


# class definitions 
# ================================================================================

class Distance(object):
    """
    A meta distance class that provides for union, intersection and difference
    opperations to compute the distance opperations of set compositions 
    """
    def __init__(self, fun=None, bbox=None):
        " set call if instantiated with a fun "
        self.fun = fun

    def __call__(self, *args, **kwds):
        " force subclasses to define this "
        return self.fun(*args, **kwds)

    def __add__(self, d):
        " return the rule for the new set "
        def dor(x,y):
            davals, dbvals = self(x,y), d(x,y)
            return where(davals < dbvals, davals, dbvals)
        return Distance(fun=dor)

    def __sub__(self, d):
        " return the distance rule for the set difference "
        def dnot(x,y):
            davals, dbvals = self(x,y), d(x,y)
            return where(davals > -dbvals, davals, -dbvals)
        return Distance(fun=dnot)

    def difference(self, d):
        " return the distance rule for the set difference. Use self - d "
        return self -  d

    def union(self, d):
        " return the rule for the new set. Use self + d"
        return self + d

    def intersect(self, d):
        " return the rule for the set intersection "
        def dand(x,y):
            davals, dbvals = self(x,y), d(x,y)
            return where(davals > dbvals, davals, dbvals)
        return Distance(fun=dand)

    def plot(self, x, y, title=None,subplot=None, levels=[-.01,.01], **kwds):
        " Make a contour plot over the given x and y "
        if subplot:
            subplot.contourf(x, y, self(x,y), levels=levels)
        else:
            plt.contourf(x, y, self(x,y), levels=levels, colors='k')
            # plt.colorbar()
            if title: plt.title(title)
        
class BoxDistance(Distance):
    " Define a box and distance function "
    def __init__(self, bbox):
        " store the bounding box as an array "
        self.bbox = array(bbox)
        self.width, self.height = bbox[1]-bbox[0], bbox[3]-bbox[2]

    def __call__(self, x, y):
        " Compute the distance outside the box "
        distance = array([self.bbox[0]-x, x - self.bbox[1],
                          self.bbox[2]-y, y - self.bbox[3]]).max(0)
        return distance

class DiskDistance(Distance):
    " define the distance to edge of a disk "
    def __init__(self, center=[0.,0.], radius=1.0):
        self.center = array(center)
        self.radius = radius
        self.bbox = array([center[0]-radius, center[0]+radius, 
                           center[1]-radius, center[1]+radius])

    def __call__(self, x, y):
        " Compute the distance outside a circle " 
        x = x - self.center[0]
        y = y - self.center[1]
        return sqrt(x**2 + y**2) - self.radius
        


# ====
# Example script
# ================================================================================

if __name__=="__main__":

    # imports

    #script
    pass

