
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
# from numpy import 
# from numexpr import 
# from tables import
import matplotlib.pyplot as plt
from enthought.mayavi import mlab 

# my imports 
from icosphere import IcoSphere, Icosahedron


# import plane_parallel


# Global variables
# ================================================================================
bgcolor = (0,0,0)
fgcolor = (1,1,1)
figsize = (1200, 900)




# make a spatial_set object
a = 0.0  
b = 10.0
nz = 100



# function and class definitions 
# ================================================================================

class SpatialSet(object):
    """ 
    Return an object with attributes for the points in space and the 
    points on the boundary. 
    """

    default_figure_number = 1

    def __init__(self, a, b, hz, ext=None):
        """
        Define points on the given interval and return. 
        
        * plan * 
        implement adaptive behavior by spacing points according to ext
        """
        # define the gridspacing
        self.hz = hz
        self.nz = (b-a) // hz
        # data
        self.points = sp.linspace(a, b, self.nz)
        self.boundary_points = sp.linspace(a, b, 2)
        self.boundary_normals = sp.array([[0., 0.,-1.0], 
                                          [0., 0., 1.0]])


        # data imbeded in three dimensions
        self.pz = self.points
        self.px = sp.zeros(self.points.shape)
        self.py = sp.zeros(self.points.shape)
        self.points_3d = sp.zeros((self.pz.size, 3))
        self.points_3d[:,2] = self.pz
        self.bz = self.boundary_points
        self.bx = sp.zeros(self.boundary_points.shape)
        self.by = sp.zeros(self.boundary_points.shape)
        self.boundary_points_3d = sp.zeros(self.boundary_normals.shape)
        self.boundary_points_3d[:,2] = self.bz


class DirectionSet(IcoSphere, Icosahedron):
    """
    Define a set of directions by specifying the approximate gritpoint spacing.
    """
    base_h = sp.sqrt(((Icosahedron.p[Icosahedron.bar[0,0]] - 
                       Icosahedron.p[Icosahedron.bar[0,1]])**2).sum())
    
    def __init__(self, h):
        " h = aproximate spacing "

        # Instantiate as an IcoSphere with appropriate resolution
        self.h = h
        self.nrefine = DirectionSet.base_h // h
        IcoSphere.__init__(self, self.nrefine) # defines self.p, self.tri, self.bar 




def get_incoming_set(spatial_set, direction_set, mincosine=1e-8):
    """
    Return an array of points (x, y, z, u, v, w) as could be used in a 
    quiver plot.
    """
    v = direction_set.p
    bp = spatial_set.boundary_points_3d
    bnorm = spatial_set.boundary_normals
    direction_cosines = sp.tensordot(bnorm, v, axes=[(-1),(-1)])
    incoming_bolean = direction_cosines < -mincosine
    
    shifts = incoming_bolean.sum(-1)
    stops = sp.cumsum(shifts)
    starts = starts = stops - shifts
    
    # initialize the incoming_set array and fill with points
    bpdim = bp.shape[-1]
    vdim = v.shape[-1]
    incoming_set_shape = (stops[-1], bpdim + vdim)
    incoming_set = sp.zeros(incoming_set_shape)
    for n, start, stop, shift in zip(sp.arange(bp.shape[0]), 
                                     starts, stops, shifts):
        incoming_set[start:stop, :bpdim] = bp[n]
        incoming_set[start:stop, bpdim:] = v[incoming_bolean[n], :]
    
    return incoming_set

def get_outgoing_set(spatial_set, direction_set, mincosine=1e-8):
    """
    Return an array of points (x, y, z, u, v, w) as could be used in a 
    quiver plot.
    """
    v = direction_set.p
    bp = spatial_set.boundary_points_3d
    bnorm = spatial_set.boundary_normals
    direction_cosines = sp.tensordot(bnorm, v, axes=[(-1),(-1)])
    incoming_bolean = direction_cosines > mincosine
    
    shifts = incoming_bolean.sum(-1)
    stops = sp.cumsum(shifts)
    starts = starts = stops - shifts
    
    # initialize the incoming_set array and fill with points
    bpdim = bp.shape[-1]
    vdim = v.shape[-1]
    incoming_set_shape = (stops[-1], bpdim + vdim)
    incoming_set = sp.zeros(incoming_set_shape)
    for n, start, stop, shift in zip(sp.arange(bp.shape[0]), 
                                     starts, stops, shifts):
        incoming_set[start:stop, :bpdim] = bp[n]
        incoming_set[start:stop, bpdim:] = v[incoming_bolean[n], :]
    
    return incoming_set

def get_interior_set(spatial_set, direction_set):
    """
    Return an array of points (x, y, z, u, v, w) as could be used in a 
    quiver plot.
    """
    v = direction_set.p
    p = spatial_set.points_3d[1:-1,:]

    # define interior set 
    interior_set_shape = (p.shape[0], v.shape[0], 6)
    interior_set = sp.zeros(interior_set_shape)
    interior_set[:, :, :3] = sp.tensordot(sp.ones(v.shape[0]), p,
                                       axes=0).transpose(1,0,2)
    interior_set[:, :, 3:] = sp.tensordot(sp.ones(p.shape[0]), v,
                                       axes=0)
    

    return interior_set




    
# for bpoint, bnormals in zip(boundary_points, boundary_normals):
            
            
            






        


        
        


# ====
# Example script
# ================================================================================

if __name__=="__main__":

    # imports

    #script


    # make a spatial_set object
    a = 0.0  
    b = 10.0
    nz = 100
    hz = .6
    hv = .25

    direction_set = DirectionSet(hv)
    spatial_set = SpatialSet(a,b,hz)
    incoming_set = get_incoming_set(spatial_set, direction_set)
    outgoing_set = get_outgoing_set(spatial_set, direction_set)
    interior_set = get_interior_set(spatial_set, direction_set)
    mlab.close(0)    
    f = mlab.figure(0, 
                    bgcolor=bgcolor, 
                    fgcolor=fgcolor, size=figsize)

    # plot the figure
    mlab.clf(0)    
    mlab.plot3d(*spatial_set.points_3d.T, tube_radius=.05)
    mlab.points3d(*spatial_set.boundary_points_3d.T, scale_factor=.3)
    mlab.quiver3d(*incoming_set.T, color=(0,0,1), scale_factor=.3)
    mlab.quiver3d(*outgoing_set.T, color=(1,0,0), scale_factor=.3)
    mlab.quiver3d(*interior_set.T, color=(0,1,0), scale_factor=.2)












