#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
import NewLap, FFT
import Plots, Constants
import matplotlib.pyplot as plt

# Return to original variable using ingterpolation
def OrigGraph(Edges,ScSoln):
    
    twopi = Constants.TWOPI
    Tsteps = np.shape(ScSoln)[0]
    Nedge = np.shape(ScSoln)[1]
    Nsamp = Constants.NSAMP
    Np1 = Nsamp + 1
    
 # Find min edgelength
    Emin = Edges[0][2]
    for j in np.arange(Nedge):
        if Edges[j][2] < Emin:
            Emin = Edges[j][2]
   
# Longest edge gets Np1 samples. 
# Other edges get Np1*l_e samples    
    
    Esamp = np.zeros(Nedge).astype(int)     
    for j in np.arange(Nedge):
        Esamp[j] = int(Nsamp*Edges[j][2])
   
    NewSoln = 1j*np.zeros((Tsteps,Nedge,Np1))
    AbSoln = np.zeros((Tsteps,Nedge,Np1))  
    print('Esamp =', Esamp )

# The stretch mapping is x to   
# phi _e = z/l_e - (1-l_e)*sin(twopi*z/l_e)/twopi 

    for t in np.arange(Tsteps):
        for k in np.arange(Nedge):
            le = Edges[k][2]
            for n in np.arange(Esamp[k]+1):
                z = n*le/Esamp[k]
                y = z/le - (1-le)*np.sin(twopi*z/le)/twopi
# m is approximate sample number in [0,1]
                m = y*Nsamp
                ind1 = int(np.floor(m))
                ind2 = int(np.ceil(m))
                wt = m - ind1
                f1 = ScSoln[t][k][ind1]
                f2 = ScSoln[t][k][ind2]
                g1 = np.abs(f1)
                g2 = np.abs(f2)
                NewSoln[t][k][n] = wt*f2 + (1-wt)*f1
                AbSoln[t][k][n] = wt*g2 +(1-wt)*g1
                if (ind1 > 0) and (wt < .5):
                    ind0 = ind1 - 1
                    f0 = ScSoln[t][k][ind0]
                    NewSoln[t][k][n] = f0*wt*(wt-1)/2 - f1*(wt-1)*(wt+1) + f2*wt*(wt+1)/2
                elif  (ind2 < np.shape(ScSoln)[2]-2) and (wt > .5):
                    wt = m - ind2
                    ind3 = ind2 + 1                        
                    f3 = ScSoln[t][k][ind3]
                    NewSoln[t][k][n] = f1*wt*(wt-1)/2 - f2*(wt-1)*(wt+1) + f3*wt*(wt+1)/2
    
    return(Esamp,NewSoln,AbSoln)