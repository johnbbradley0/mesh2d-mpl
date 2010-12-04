
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
from matplotlib.delaunay import delaunay
from enthought.mayavi import mlab
from scipy.spatial.distance import pdist

# from numpy import 
# from numexpr import 
# from tables import

# my imports


# Global variables
# ================================================================================




    


# function and class definitions 
# ================================================================================

def get_points():
    " Define the 12 verticies on the z-axis aligned icosahedron "

    # Define the verticies with the goldon ratio
    a = (1 + sp.sqrt(5))/2.0   # goldon ratio
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

    return p[reorder_index, :]


def get_barymat(n):
    """
    Define the matrix that will refine points on a triangle 
    
    """
    numrows = n*(n+1)/2

    # define the values that will be needed 
    ns = sp.arange(n)
    vals = ns / float(n-1)
 
    # initialize array
    bcmat = sp.zeros((numrows, 3))#sp.arange(numrows*3).reshape((numrows, 3))

    # loop over blocks to fill in the matrix
    shifts = sp.arange(n,0,-1)
    starts = sp.zeros(n, dtype=int); 
    starts[1:] += sp.cumsum(shifts[:-1]) # starts are the cumulative shifts
    stops = starts + shifts              
    for n_, start, stop, shift in zip(ns, starts, stops, shifts):
        bcmat[start:stop, 0] = vals[shift-1::-1]
        bcmat[start:stop, 1] = vals[:shift]
        bcmat[start:stop, 2] = vals[n_]

    return bcmat


    



    


class Icosahedron(object):
    """ 
    The verticies of an icosahedron, together with triangles, edges, and 
    triangle midpoints and edge midpoints.  The class data stores the 
    """
    
    # define points (verticies)
    p = get_points()
    px, py, pz = p[:,0], p[:,1], p[:,2]
    
    # define triangles (faces)
    tri = sp.array([[1,2,3,4,5,6,2,7,2,8,3, 9,10,10,6, 6, 7, 8, 9,10],
                    [2,3,4,5,1,7,1,8,8,9,9,10, 5, 6,1,11,11,11,11,11],
                    [0,0,0,0,0,1,7,2,3,3,4, 4, 4, 5,5, 7, 8, 9,10, 6]]).transpose()
    trimids = (p[tri[:,0]] + p[tri[:,1]] + p[tri[:,2]]) / 3.0    

    # define bars (edges)
    bar = list()
    for t in tri:
        bar += [sp.array([i,j]) for i, j in [t[0:2], t[1:], t[[2, 0]]] if j>i]
    bar = sp.array(bar)
    barmids = (p[bar[:,0]] + p[bar[:,1]]) / 2.0
    
    


# ====
# Example script
# ================================================================================
def triangulate_bary(bary):
    """
    triangulate a barycentric triangle using matplotlib.
    return (bars, triangles)
    """
    x, y = sp.cos(-sp.pi/4.)*bary[:,0] + sp.sin(-sp.pi/4.)*bary[:,1] , bary[:,2]
    dely = delaunay(x,y)
    return dely[1], dely[2]


def get_triangulation(n, ico=Icosahedron()):
    """
    Compute the triangulation of the sphere by refineing each face of the 
    icosahedron to an nth order barycentric triangle.  There are two key issues
    that this routine addresses.
    
    1) calculate the triangles (unique by construction)
    2) remove non-unique nodes and edges

    """
    
    verts = sp.array([ico.p[ico.tri[:,0]],
                      ico.p[ico.tri[:,1]],
                      ico.p[ico.tri[:,2]]])

    bary = get_barymat(n)
    newverts = sp.tensordot(verts, bary,  axes=[(0,), (-1,)]).transpose(0,2,1)    
    numverts = newverts.shape[1]
    if newverts.size/3 > 1e6: print("newverts.size/3 is high: {0}".format(
            newverts.size/3))
    flat_coordinates = sp.arange(newverts.size/3).reshape(20, numverts)


 
    barbary, tribary = triangulate_bary(bary)
    newtri = sp.zeros((20, tribary.shape[0], 3), dtype=int)
    newbar = sp.zeros((20, barbary.shape[0], 2), dtype=int)
    for i in range(20):
        for j in range(3):
            newtri[i, :, j] = flat_coordinates[i, tribary[:,j]]
            if j < 2: newbar[i, :, j] = flat_coordinates[i, barbary[:,j]]
            


    newverts = newverts.reshape(newverts.size/3, 3)
    newtri = newtri.reshape(newtri.size/3, 3)
    newbar = newbar.reshape(newbar.size/2, 2)

    # normalize verticies
    scalars = sp.sqrt((newverts**2).sum(-1))
    newverts = (newverts.T / scalars).T


    # remove repeated verticies 
    aux, iunique, irepeat = sp.unique(sp.dot(newverts//1e-8, 100*sp.arange(1,4,1)), 
                                      return_index=True, return_inverse=True)

    univerts = newverts[iunique]
    unitri = irepeat[newtri]
    unibar = irepeat[newbar]
    mid = .5 * (univerts[unibar[:,0]] + univerts[unibar[:,1]])
    aux, iu  = sp.unique(sp.dot(mid//1e-8, 100*sp.arange(1,4,1)), return_index=True) 
    unimid = mid[iu]
    unibar = unibar[iu,:]

    return univerts, unitri, unibar


class IcoSphere(Icosahedron):
    """
    """
    def __init__(self, n):
        """
        define an icosahedron based discritization of the sphere
        n is the order of barycentric triangles used to refine each 
        face of the icosaheral base mesh.
        """
        self.p, self.tri, self.bar = get_triangulation(n, Icosahedron)
        
        


# newtri = sp.zeros((20,tribary.shape[0], 3))
 
# for i in range(20):
#     newtri[i,:,:] = tribary
    





























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
