# Create all points on a lattice discretisation of a hypercube of N dimensions

import numpy as np
import itertools

# define hypercube
N = 3
splitFiles = 1
#latticePoints = np.linspace(0.5,2,2)
latticePoints = [0.1,0.5,1.0,1.5,2.0,2.5,3.0]
#latticePoints[0] += 0.1  # maybe not really zero, just put it small enough
allDimSameSpacing = True

if allDimSameSpacing:
    toIter = [latticePoints for i in xrange(N)]
    wholeLatticeWOIdx = np.array(list(itertools.product(*toIter)))
else:
    # Define later, if needed.
    pass
'''
# insert index column as the first column
wholeLattice = np.zeros((wholeLatticeWOIdx.shape[0],wholeLatticeWOIdx.shape[1]+1))
wholeLattice[:,0] = np.arange(wholeLatticeWOIdx.shape[0])
wholeLattice[:,1:] = wholeLatticeWOIdx[:,:]


if len(wholeLattice)%splitFiles:
    print "WARNING: files not equally split or splitting may have gone wrong."

n = len(wholeLattice)/splitFiles
for i in range(splitFiles):
    np.savetxt("conductanceDataTest%d.txt"%i,wholeLattice[i*n:(i+1)*n])
'''

thefile = open('conductanceDataTest.txt', 'w')
for i in xrange(len(wholeLatticeWOIdx)):
    thefile.write("%d"%i)
    for a in wholeLatticeWOIdx[i]:
        thefile.write("\t%f"%a)
    thefile.write("\n")


