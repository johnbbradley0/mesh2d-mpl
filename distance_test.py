
#
# Author: William G.K. Martin (wgm2111@cu where cu=columbia.edu)
# copyright (c) 2010
# liscence: BSD style
#
# ================================================================================


# standard python imports
from __future__ import print_function

# third party imports 
from scipy import 
from numpy import 
from numexpr import 
from tables import
import nose.tools as nt
import numpy.testing as npt

# my imports 
import distance as dist

# Global variables
# ================================================================================

# making three test circles 
disk0 = dist.DiskDistance([.5, sqrt(3)/2], 1.25)
disk1 = dist.DiskDistance([.5, -sqrt(3)/2], 1.25)
disk2 = dist.DiskDistance([-1, 0], 1.25)

in0out12 = []


# function and class definitions 
# ================================================================================



# ====
# Example script
# ================================================================================

if __name__=="__main__":

    # imports

    #script
    pass

