
"""
"""
import numpy as np
import Firststeps, Constants,Sigmas 
import Askar 
import NewLap
import CombLap, QGvec
import PiCase, Util, FFT
import InitData, PIFT, Plots
import matplotlib.pyplot as plt
import Map, Adjust
import time

Edges,Adjac0,Eorder,Edir = Firststeps.firststuff()

Graphcase = 22

# Several functions result from the conversion of all edges to length one
sigrt,sigsq,sigpot = Sigmas.Sigdef(Edges)

# Find the adjacency matrix, eigenvalues and eigenvectors 
Adj,eigenvalues,eigenvectors = CombLap.CombEval(Edges)
print('Combinatorial eigenvalues')
print(eigenvalues)

# Make a list showing each vertex, its adjacent vertices, and the corresponding
# edge number.
Vlist = Util.XVertlist(Edges)

# Get quantum graph fundamental frequencies and coefficients
# for quantum graph eigenfunctions coming from the combinatorial eigenvalues. 
QGFrqs,QGEvec = QGvec.QGvec1(Edges,eigenvalues,eigenvectors)

# Get 'Pi' eigenvectors', which were exceptional cases not handled in the previoius step.
PiFrqs,PiEvecs = PiCase.PiComps(Adj,Edges)

Combnum = np.shape(QGFrqs)[0]
Pinum = np.shape(PiFrqs)[0]
print('Pinum = ',Pinum,'Combnum =',Combnum, 'Total = ',Pinum+Combnum)

# Merge two eigenvalue lists
Frqs,Evecs = Util.Merge(QGFrqs,QGEvec,PiFrqs,PiEvecs)

print('Merged omegas = sqrt(lambda)')
print(Frqs)

Nedge = np.shape(Evecs)[0]
Nfrqs = np.shape(Frqs)[0]

Finit = InitData.SchrdData(sigrt,Edges,Nedge,Frqs,Evecs)

# Solve equation 
Tfinal = float(input('Enter final time '))
Tfinal = Tfinal*Constants.TIMESCALE

Tsteps = int(input('Select number of time steps '))
thyme = Tfinal/Tsteps

Nsamp = Constants.NSAMP
Soln = 1j*np.zeros((Tsteps+1,Nedge,Nsamp+1))
Soln[0][:][:] = Finit 

for stepnum in np.arange(1,Tsteps+1):
    Soln = Askar.E1(Frqs,Evecs,stepnum,thyme,sigsq,sigpot,Soln)
    if np.mod(stepnum,50) == 1:
        print('Completed step ',stepnum)

# Rescale solution 
ScSoln = Map.SMap(sigrt,Soln)

# Return to orignal edge variables.
Esamp,NewSoln,AbSoln = Adjust.OrigGraph(Edges,ScSoln)


# Plot absolute value squared of solution
Plots.Loopy(Eorder,Edir,Esamp,AbSoln)


print('Finished')
