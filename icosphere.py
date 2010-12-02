
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
import scipy.linalg as la 

# from numpy import 
# from numexpr import 
# from tables import

# my imports


# Global variables
# ================================================================================

# a = (1 + sp.sqrt(5))/2.0
# p = sp.array([[a,-a,-a, a, 1, 1,-1,-1, 0, 0, 0, 0],
#            [0, 0, 0, 0, a,-a,-a, a, 1, 1,-1,-1],
#            [1, 1,-1,-1, 0, 0, 0, 0, a,-a,-a, a]]).transpose()
# pnorms = sp.sqrt((p**2).sum(1))
# p = p / pnorms[0]
# ang = sp.arctan(p[0,0] / p[0,2])
# ca, sa = sp.cos(ang), sp.sin(ang)
# rotation = sp.array([[ca, 0, -sa], [0, 1.0, 0], [sa, 0, ca]])
# p = sp.inner(rotation, p).transpose()


# weights = sp.arange(p.shape[0], dtype=float)

# # reorder the points by finding the tip, upper, lower, bot
# reorder_index = [0, 3, 4, 8, -1, 5,-2, -3, 7, 1, 6, 2]
# p = p[reorder_index, :]

# # define triangles
# tri = sp.array([[1,2,3,4,5,6,2,7,2,8,3, 9,10,10,6, 6, 7, 8, 9,10],
#                 [2,3,4,5,1,7,1,8,8,9,9,10, 5, 6,1,11,11,11,11,11],
#                 [0,0,0,0,0,1,7,2,3,3,4, 4, 4, 5,5, 7, 8, 9,10, 6]]).transpose()
# trimids = (p[tri[:,0]] + p[tri[:,1]] + p[tri[:,2]]) / 3.0
# triweights = sp.arange(tri.shape[0], dtype=float)


# # define bars (edges)
# bars = list()
# for t in tri:
#     bars += [sp.array([i,j]) for i, j in [t[0:2], t[1:], t[[2, 0]]] if j>i]
# bars = sp.array(bars)
# barmids = (p[bars[:,0]] + p[bars[:,1]]) / 2.0

class Icosahedron(object):
    """ 
    The verticies of an icosahedron, together with triangles, edges, and 
    triangle midpoints and edge midpoints 
    """
    def __init__(self):
        # define the points (verticies)
        a = (1 + sp.sqrt(5))/2.0
        p = sp.array([[a,-a,-a, a, 1, 1,-1,-1, 0, 0, 0, 0],
                      [0, 0, 0, 0, a,-a,-a, a, 1, 1,-1,-1],
                      [1, 1,-1,-1, 0, 0, 0, 0, a,-a,-a, a]]).transpose()
        p = p / sp.sqrt((p**2).sum(1))[0]
        # rotate top point to the z-axis
        ang = sp.arctan(p[0,0] / p[0,2])
        ca, sa = sp.cos(ang), sp.sin(ang)
        rotation = sp.array([[ca, 0, -sa], [0, 1.0, 0], [sa, 0, ca]])
        p = sp.inner(rotation, p).transpose()    
        # reorder in a downward spiral
        reorder_index = [0, 3, 4, 8, -1, 5,-2, -3, 7, 1, 6, 2]
        p = p[reorder_index, :]
    
        # define triangles
        tri = sp.array([[1,2,3,4,5,6,2,7,2,8,3, 9,10,10,6, 6, 7, 8, 9,10],
                        [2,3,4,5,1,7,1,8,8,9,9,10, 5, 6,1,11,11,11,11,11],
                        [0,0,0,0,0,1,7,2,3,3,4, 4, 4, 5,5, 7, 8, 9,10, 6]]).transpose()
        trimids = (p[tri[:,0]] + p[tri[:,1]] + p[tri[:,2]]) / 3.0
        triweights = sp.arange(tri.shape[0], dtype=float)
        
        # define bars (edges)
        bar = list()
        for t in tri:
            bar += [sp.array([i,j]) for i, j in [t[0:2], t[1:], t[[2, 0]]] if j>i]
        bar = sp.array(bar)
        barmids = (p[bar[:,0]] + p[bar[:,1]]) / 2.0
    
        # store data
        self.p = p
        self.px, self.py, self.pz = p[:,0], p[:,1], p[:,2]
        self.tri, self.trimids = tri, trimids 
        self.bar, self.barmids = bar, barmids
    

    


# function and class definitions 
# ================================================================================



# ====
# Example script
# ================================================================================

if __name__=="__main__":

    # imports
    from enthought.mayavi import mlab
    
    #script

    # define icosahedron object
    ico = Icosahedron()

    # define weights for showing index order
    weights = sp.arange(ico.px.size, dtype=float)
    triweights = sp.arange(ico.tri.shape[0], dtype=float)

    # make a plot of the icosahedron
    mlab.close(0)
    mlab.figure(0,bgcolor=(.1, .1, .1), fgcolor=(1., 1., 1.), size=(800, 700))
    mlab.plot3d(ico.px, ico.py, ico.pz, weights)
    mlab.points3d(ico.px, ico.py, ico.pz, weights, scale_mode='none')
    mlab.triangular_mesh(ico.px, ico.py, ico.pz, ico.tri, scalars=weights)
    mlab.points3d(*ico.barmids.T, scale_factor=.07)
    mlab.points3d(*ico.trimids.T, scale_factor=.1)
    mlab.colorbar()
