# Introduction #

The meshing routines here are being designed to create and store adaptive unstructured meshes for use in radiative transfer calculations.  The project is a one man hack job that uses the well known python libraries scipy, numpy, and matplotlib.


# Details #

Implementing the mesh will be fairly strait forward and accomplished by the following procedure:

## define uniform points in the domain's bounding box ##
  1. remove outside points with 'd'
  1. remove low probability points with 'h'
  1. add fixed points specified by the user

## define a grid with the points object ##
  1. run delaunay for getting triangles
  1. extract bars and centers from the triangulation
  1. plot the grid

## iterate over force based adjustments ##
  1. re-triangulate if the points moved enough
  1. move points with push-only force vector
    * zero force on fixed points
  1. project out-of-bounds points to the boundary
## return a mesh object with access to the following data ##
  1. points in 2D
  1. bars
    * bar\_normals
    * bar\_lengths
  1. triangles
  1. triangle\_centers
  1. triangle\_area
  1. plot\_mesh (points, bars and triangles)


This procedure sketches the instantiation of a mesh and will rely on user provided data and function definitions.

# input: #
  1. d2b(points) = "signed distance to the boundary"
  1. h\_local(points) = "local grid spacing"
  1. pfix = "grid points that are fixed"
  1. h0 = "starting resolution"
  1. bbox = [xmin, xmax, ymin, ymax]

# constants: #
  1. dptol = "minimum interesting step length"
  1. Fscale = "scale desired length by 1.2 for force calculation"
  1. dt = "time step"
  1. geps = "triangle center remove distance, .001\*h0"
  1. deps = "length for numerical differentiation" sqrt(eps)