
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================

"""

Intended to by used to create a mesh object for radiative transfer simulations 
in the plane parallel geometry.

The only Object that should be needed for external purposes is "Domain".  Create
a domain instance as follows.

**example**
>>> a, b, hz, hv = 0.0, 10.0, .6, .3
>>> domain_set = Domain(a, b, hz, hv)

** see also **
domain.quiver3d and domain.plot_function

"""


# standard python imports
from __future__ import print_function


# third party imports
import scipy as sp
import matplotlib.pyplot as plt
from enthought.mayavi import mlab 

# my imports 
from icosphere import IcoSphere, Icosahedron

# my imports 


# Global variables
# ================================================================================



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
        self.icosphere = IcoSphere(self.nrefine) # defines self.p, self.tri, self.bar 
        self.directions = self.icosphere.p
        

# functions for computing the point_sets for different submanifolds
def get_incoming_set(spatial_set, direction_set, mincosine=1e-8):
    """
    Return an array of points (x, y, z, u, v, w) as could be used in a 
    quiver plot.
    """
    v = direction_set.directions
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
    v = direction_set.directions
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

def get_domain_set(spatial_set, direction_set):
    """
    Return an array of points (x, y, z, u, v, w) as could be used in a 
    quiver plot.
    """
    v = direction_set.directions
    p = spatial_set.points_3d#spatial_set.points_3d[1:-1,:]

    # define interior set 
    interior_set_shape = (p.shape[0], v.shape[0], 6)
    interior_set = sp.zeros(interior_set_shape)
    interior_set[:, :, :3] = sp.tensordot(sp.ones(v.shape[0]), p,
                                       axes=0).transpose(1,0,2)
    interior_set[:, :, 3:] = sp.tensordot(sp.ones(p.shape[0]), v,
                                       axes=0)
    newshape = (interior_set.size // 6, 6)
    # newshape = (interior_set_shape[0]*interior_set_shape[1], 6)
    return interior_set.reshape(newshape)


# 
# Define domain

class Domain(object):
    """
    Store the spatial points, directions, and three inportant point_sets.
    
    incoming_set - point, direction
    outgoing_set
    domain_set

    """
    def __init__(self, a, b, hz, hv):
        """ 
        Compute all the pointsets needed for a plane-parallel radiative transfer 
        domain. 
        """
        # instantiate direction and spatial sets
        self.spatial_set = SpatialSet(a,b,hz)
        self.direction_set = DirectionSet(hv)
        
        # store the spatial points and directions  
        self.points = self.spatial_set.points
        self.boundary_points = self.spatial_set.boundary_points
        self.points_3d = self.spatial_set.points_3d
        self.boundary_points_3d = self.spatial_set.boundary_points_3d
        self.directions = self.direction_set.directions
        
        # store the point_sets in 3d
        self.domain_set_3d = get_domain_set(self.spatial_set, self.direction_set)
        self.incoming_set_3d = get_incoming_set(self.spatial_set, self.direction_set)
        self.outgoing_set_3d = get_outgoing_set(self.spatial_set, self.direction_set)
        
        # store views of the 3d arrays as point_sets in 1D
        for name in ["domain_set", "incoming_set", "outgoing_set"]:
            name3d = name+"_3d"
            setattr(self, name, getattr(self, name3d)[..., 2:]) 
            

    def plot_function(self, fun, **kwds):
        """
        Plot a function of function values. 
        """
        fun_values = fun(self.points)
        args = tuple(self.points_3d.T) + (fun_values,)
        mlab.plot3d(*args, **kwds)

        fun_boundary_values = fun(self.boundary_points)
        args = tuple(self.boundary_points_3d.T) + (fun_boundary_values,)
        mlab.points3d(*args, scale_factor=.2, scale_mode='none')



    def quiver3d(self, name, **kwds):
        " Plot the domain with a quiver plot "
        names = ['domain_set_3d', 'incoming_set_3d', 'outgoing_set_3d'] 
        if name in names:
            mlab.quiver3d(*getattr(self, name).T, **kwds)
        else:
            print("WARNING: name = {0} not in {1}".format(
                    name, str(names)))


# ====
# Example script
# ================================================================================

if __name__=="__main__":        
    pass
