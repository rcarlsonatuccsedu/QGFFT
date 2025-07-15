#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 08:42:52 2023

@author: robertcarlson

"""

import numpy as np

# Input: The edge sample number and FFT size parameter,
# the combinatorial eigenvalues and eigenvectors
# and the n^2pi^2 eigenvalues and eigenvectors
# Output: A merged eigenvalue list and corresponding eigenvector list
def Merge(QGEvals,QGEvecs,PiEvals,PiEvecs):

    import matplotlib.pyplot as plt
    import Constants
     
    Nsamp = Constants.NSAMP
    QGnum = np.shape(QGEvals)[0]
    Pinum = np.shape(PiEvals)[0]
    
    Nedges = np.shape(QGEvecs)[0]
    
    Evalcnt = QGnum+Pinum
    pntr1 = 0
    pntr2 = 0
    
    Frqs = np.zeros(Evalcnt)
    Evecs = 1j*np.zeros((Nedges,Evalcnt,2))
    
    case = 0
    for k in np.arange(Evalcnt):
        if pntr2 == Pinum: case = 1
        if pntr1 == QGnum: case = 2
        
        if case == 0:
            if QGEvals[pntr1] < PiEvals[pntr2]: 
                case = 1
            else: case = 2
            
        if case == 1:
            Frqs[k] = QGEvals[pntr1]
            for j in np.arange(Nedges):
                Evecs[j][k][0] = QGEvecs[j][pntr1][0]
                Evecs[j][k][1] = QGEvecs[j][pntr1][1]
            pntr1 += 1    
        if case == 2:
            Frqs[k] = PiEvals[pntr2]
            for j in np.arange(Nedges):
                Evecs[j][k][0] = PiEvecs[j][pntr2][0]
                Evecs[j][k][1] = PiEvecs[j][pntr2][1]
            pntr2 +=1    
        case = 0 
      
    # Normalize Evecs 
    FF = np.sqrt(2/Nsamp)
    for k in np.arange(Evalcnt):
        for j in np.arange(Nedges):
            Evecs[j][k][0] = FF*Evecs[j][k][0]
            Evecs[j][k][1] = FF*Evecs[j][k][1]
      
    return (Frqs,Evecs)

# Input: Edge list
# Output: Prints a list of adjacent vertex pairs and the associated 
# edge numbers
def XVertlist(Edges):
    
    Nedge = np.shape(Edges)[0]
    Vlist = np.zeros((2*Nedge,3))
 
    vcnt = 0
    
    for j in np.arange(Nedge):
        if Edges[j][0] > vcnt:
            vcnt = Edges[j][0]
        if Edges[j][1] > vcnt:
            vcnt = Edges[j][1]
            
    pntr = 0        
    for j in np.arange(vcnt+1):
            for k in np.arange(Nedge):
                if Edges[k][0] == j:
                    Vlist[pntr][0] = j
                    Vlist[pntr][1] = Edges[k][1]
                    Vlist[pntr][2] = k
                    pntr += 1
                if Edges[k][1] == j:
                    Vlist[pntr][0] = j
                    Vlist[pntr][1] = Edges[k][0]
                    Vlist[pntr][2] = k
                    pntr += 1
                
    for j in np.arange(2*Nedge):
        print('v1 = ',Vlist[j][0],'v2 = ',Vlist[j][1], 'edge = ',Vlist[j][2])



    return(Vlist)