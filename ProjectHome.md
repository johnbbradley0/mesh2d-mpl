To enable adaptivity in 2D and 3D radiative transfer simulations, this mesh code is based on a meshing algorithms given by persson\_simple\_2002 and is under heavy [development](http://code.google.com/p/mesh2d-mpl/wiki/development).  The code will be designed to handle geometric and solution based adaptivity through the use of two functions: one computing the distance from a point to the boundary and a second that gives the desired length resolution at each point.  Together they are used to adapt the grid resolution.