#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 08:42:52 2023

@author: robertcarlson
"""
import numpy as np
import PiFunks
# This function generatest the quantum graph frequency and eigenfunction data
# associated to combinatorial eigenvalues (n\pi )^2.

def PiComps(Adj, Edges):

    Nedge = np.shape(Edges)[0]
 
    freq = np.pi
    Eigvecspi = PiFunks.PiEqns(freq, Adj, Edges)
    Veccnt1 = np.shape(Eigvecspi)[1]
    freq = 2*np.pi
    Eigvecs2pi = PiFunks.PiEqns(freq, Adj, Edges)
    Veccnt2 = np.shape(Eigvecs2pi)[1]
 
    Veccnt0 = 1

    Veccnt = Veccnt0 + Veccnt1 + Veccnt2

    PiFrqs = np.zeros(Veccnt)
    PiEvecs0 = 1j*np.zeros((Nedge, Veccnt, 2))
    PiEvecs = 1j*np.zeros((Nedge, Veccnt, 2))

# The first eigenvalue is zero, with constant eigenfunction
# Normalize the constant eigenfunction

    PiFrqs[0] = 0
    Nvert = int(np.shape(Adj)[0])   
    Degsum = 0
    for j in np.arange(Nvert):
        for k in np.arange(Nvert):
             Degsum = Degsum + Adj[j][k]
    normval = 1/np.sqrt(Degsum)
    for k in np.arange(Nedge):
        PiEvecs0[k][0][0] = normval
        PiEvecs0[k][0][1] = 0

# Next include the eigenvalues pi and 2*pi according to their
# multiplicties     
    pntr = 1
    for j in np.arange(Veccnt1):
        PiFrqs[pntr] = np.pi
        for k in np.arange(Nedge):
            PiEvecs0[k][pntr][0] = Eigvecspi[k][j][0]
            PiEvecs0[k][pntr][1] = Eigvecspi[k][j][1]
        pntr += 1

    for j in np.arange(Veccnt2):
        PiFrqs[pntr] = 2*np.pi
        for k in np.arange(Nedge):
           PiEvecs0[k][pntr][0] = Eigvecs2pi[k][j][0]
           PiEvecs0[k][pntr][1] = Eigvecs2pi[k][j][1]
        pntr += 1

# Change to complex exponential coefficients
    for m in np.arange(Veccnt):
        for n in np.arange(Nedge):
            PiEvecs[n][m][0] = 0.5*(PiEvecs0[n][m][0] - 1j*PiEvecs0[n][m][1])
            PiEvecs[n][m][1] = 0.5*(PiEvecs0[n][m][0] + 1j*PiEvecs0[n][m][1])      
          
    return(PiFrqs,PiEvecs)

