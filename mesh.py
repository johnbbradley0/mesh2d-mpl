
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
from scipy import linspace, array, meshgrid, nonzero, sqrt, cross, exp, rand, ones
from scipy import where, sparse, sin, pi
import matplotlib.delaunay as delaunay
from numpy import finfo, float64

# et_dely = delaunay.delaunay
import matplotlib.pyplot as plt
from numexpr import evaluate
# from tables import

# my imports 
from distance_example import dtotal, dcloud


# Global variables
# ================================================================================
d = dtotal-dcloud
bound = [-2.1,2.1,-.1,2.1]
xn = 500; 
yn = 500; 


h0 = .03251
geps = h0 * .001 #min(bound[1]-bound[0], bound[3]-bound[2])
latmosphere = lambda x,y: h0 * ones(x.shape) * (
    (1.0 + 0.3*sin(2*pi*y) + .3*sin(4*pi*x)))#* exp(.55*y) 
lcloud = lambda x,y: .95 * h0 *ones(x.shape)

# function and class definitions 
# ================================================================================



# delaunay triangulation 
def get_dely(d, nodes):
    """
    get the delaunay tesselation
    return centers, triangles, neighbors, bars, barmids
    """
    x,y = nodes[0].squeeze(), nodes[1].squeeze()
    centers, bars, tris, neighs = delaunay.delaunay(*nodes)#x,y)#*nodes)

    # remove outside triangles
    kin = nonzero( d(*centers.T)<geps)
    centers = centers[kin].T
    tris = tris[kin]
    neighs = neighs[kin]
    
    # remove outside bars
    barmids = .5 * (nodes[:,bars[:,0]] + nodes[:,bars[:,1]])
    # kin = nonzero(d(*barmids)<geps)
    bars = bars[d(*barmids)<geps]

    barmids = .5 * (nodes[:,bars[:,0]] + nodes[:,bars[:,1]])
    return centers, tris, neighs, bars, barmids

# node rejection functions
nodes_reject_oob = lambda d,nodes : nodes[:, d(*nodes)<geps]

def get_domain_area(nodes, tris):
    " calcualte the area in the triangles "
    u = nodes[:, tris[:,1]] - nodes[:,tris[:,0]]
    v = nodes[:, tris[:,2]] - nodes[:,tris[:,0]]
    area = 1.0 / 2.0 * cross(u,v,axis=0).sum()
    return area

def select_nodes(lfun, nodes, dely):
    " Get the right number of nodes for l_fun "
    # names
    centers, triangles = dely[0:2]
    # define the probability of keeping each node
    nstart = nodes.shape[1]
    area = get_domain_area(nodes, dely[1])
    lvals = lfun(*nodes)**2
    nend = 2.0 / sqrt(3) * area / lvals.mean()
    lvalsinv = 1.0 / lvals
    pkeep = nend * lvalsinv / (lvalsinv).sum() 

    # keep only the nodes that are
    return nodes[:, pkeep > rand(pkeep.size)]
   
# define a routine to move nodes away from eachother
def move_nodes(d, lfun, nodes):#, dely):
    " Move nodes one timestep (projecting lost points to the boundary) "

    # get delauney 
    dely = get_dely(d, nodes)

    # force constant 
    dt = .1520015
    deps = h0 * sqrt(finfo(float64).eps)
    restlength_factor = 1.400025
    
    # bars and midpoints
    bars, barmids = dely[-2:]
    barvecs = nodes[:,bars[:,1]] - nodes[:,bars[:,0]]
    barls = sqrt((barvecs**2).sum(0))
    u = barvecs / barls       # unit vectors
    
    # force from each bar
    restlen = restlength_factor * lfun(*barmids)

    # print('restlen: {0} \nbarls : {1}'.format(restlen.shape, barls.shape))
    logic = restlen > barls
    f = where(logic, restlen-barls, 0.0)
    # f = where(f<h0/2.0, f, h0/2.0)
    # 
    ns = nodes.shape
    spmat = sparse.csc_matrix
    # print(ns)
    # print(u[0].shape)
    # print(f.shape)
    # print(bars[:,0].shape)

    dp =  (-spmat((u[0]*f, (0*bars[:,0], bars[:,0])), shape=ns).todense() +
            spmat((u[0]*f, (0*bars[:,1], bars[:,1])), shape=ns).todense())
    dp += (-spmat((u[1]*f, (0*bars[:,0]+1, bars[:,0])), shape=ns).todense() +
            spmat((u[1]*f, (0*bars[:,1]+1, bars[:,1])), shape=ns).todense())

    nodes = array(nodes + dt*dp)

    # project boundary points back into the domain
    d_ = d(*nodes)
    ix = nonzero(d_>0)
    some_out = True
    count = 0
    while some_out:
        gradx = 1.0/deps * (d(nodes[0,ix]+deps, nodes[1,ix]) - d_[ix])
        grady = 1.0/deps * (d(nodes[0,ix], nodes[1,ix] + deps) - d_[ix])
        norm = sqrt(gradx**2 + grady**2)
        nodes[0,ix] -= d_[ix]*gradx / norm
        nodes[1,ix] -= d_[ix]*grady / norm
        d_ = d(*nodes)
        ix = nonzero(d_>geps)
        some_out = ix[0].size
        count+=1
        if count>5: #raise ValueError("counted "+str(ix[0].size)+" nodes oob")
            print("counted ",str(ix[0].size)," nodes oob")
            break
    return nodes

    


# ====
# Example script
# ================================================================================


# make surounding grid
x = linspace(bound[0],bound[1],xn)
y = linspace(bound[2],bound[3],yn)
x,y = meshgrid(x,y); y[:, ::2] = y[:, ::2] + .5 * (y[1,0] - y[0,0]) 
plt.close()    
# plt.figure(figsize=(13,17),facecolor='w')

for i, d, lfun in zip(range(10), 
                      [dtotal, dcloud, dtotal-dcloud],
                      [latmosphere, lcloud, latmosphere]):

    # make starting nodes
    nodes = array([x.flatten(), y.flatten()])
    nodes = nodes_reject_oob(d, nodes)
    Nnodes = nodes.shape[1]
    
    # make a delaunay
    dely = get_dely(d, nodes)

    # reject the nodes
    nodes = select_nodes(lfun, nodes, dely)
    dely = get_dely(d, nodes)
    
    # move nodes once
    for j in range(200):
        nodes = move_nodes(d,lfun, nodes) # dely)
    dely = get_dely(d, nodes)
    # plt.subplot(3,1,i+1)
    # d.plot(x,y)
    # plt.triplot(nodes[0],nodes[1], dely[1])
    if i ==1: 
        cloudnodes = nodes
        clouddely = dely
    elif i==2:
        clearnodes = nodes
        cleardely = dely

# show the plot
# plt.show()

# figure showing both meshes on the same plot
plt.close(10)
plt.figure(10, figsize=(25,10), facecolor='w')
(dtotal-dcloud).plot(x,y)
z = lcloud(*clouddely[0])
plt.triplot(cloudnodes[0], cloudnodes[1], clouddely[1]); plt.colorbar()
plt.triplot(clearnodes[0], clearnodes[1], cleardely[1])
# plt.tricontourf(cloudnodes[0], cloudnodes[1], clouddely[1], lcloud(*cloudnodes))
plt.tricontourf(clearnodes[0], clearnodes[1], cleardely[1], latmosphere(*clearnodes))
plt.show()


if __name__=="__main__":

    # imports

    #script
    pass

