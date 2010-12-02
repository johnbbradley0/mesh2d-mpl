
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

# from numexpr import 
# from tables import
# from matplotlib import


# my imports 
from icosphere import Icosahedron 

# Global variables
# ================================================================================


# define icosahedron object
ico = Icosahedron()

# define weights for showing index order
weights = sp.arange(ico.px.size, dtype=float)
triweights = sp.arange(ico.tri.shape[0], dtype=float)

# function and class definitions 
# ================================================================================



# ====
# Output 
# ================================================================================

if __name__=="__main__":

    # imports
    from enthought.mayavi import mlab

    # make a plot of the icosahedron
    mlab.figure(0,bgcolor=(.1, .1, .1), fgcolor=(1., 1., 1.), size=(800, 700))
    mlab.clf()
    mlab.plot3d(ico.px, ico.py, ico.pz, weights, tube_radius = .05)
    mlab.points3d(ico.px, ico.py, ico.pz, weights, 
                  scale_mode="none", scale_factor=0.3)
    mlab.triangular_mesh(ico.px, ico.py, ico.pz, ico.tri, scalars=weights)
    mlab.points3d(*ico.barmids.T, color=(0.25,0.25,0.25), scale_factor=.07)
    mlab.points3d(*ico.trimids.T, color=(0.25,0.25,0.25), scale_factor=.1)
    mlab.colorbar()


