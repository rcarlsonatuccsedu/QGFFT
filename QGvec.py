#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 15:23:13 2024

@author: robertcarlson
"""

# The equilateral quantum graph eigenfunctions and fundamental frequencies
# are computed

import numpy as np

def QGvec1(Edges,eigenvalues,eigenvectors):
    
# For each combinatorial eigenvalue nu there are qg eigenvalues lambda
# with nu = 1 - cos(sqrt(lambda))
# Exceptions are lambda = 0,pi,2*pi, which correspond to nu = 0,2   
# Values 0 <= sqrt(lambda ) <= 2pi are fundamental frequencies

# Define a tolerance for identifying special numbers nu = 0,2  
    TOL = 1.0e-14
 
# Get the number of combinatorial eigenvalues and eigenvectors
    Ncomb = np.size(eigenvalues)

# Record the qg fundamental frequncies coming from combinatorial eigenvalues
# Make an array containing the combinatorial eigenvalues and 
# explicitly recording their original index for access to combinatorial
# eigenvectors.
# The eigenvalue zero is treated later.
    qgcnt = 0
    for k in np.arange(Ncomb):
        x = eigenvalues[k]      
        if (abs(x) > TOL) and abs(2.0 - x) > TOL: 
            qgcnt += 2

    qgfrq = np.zeros((qgcnt,2))
 
# Combinatorial eigenvalues come in increasing order.
# Each one generally contributes two qg fundamental frequencies
# between zero and 2pi.  Other than zero, the exception is when
# eigenvalues[k] has value 2.  These correspond to frequencies
# equal to pi or two pi.

    valcnt = 0
    for k in np.arange(Ncomb):
        x = eigenvalues[k]      
        if (abs(x) > TOL) and abs(2.0 - x) > TOL: 
            qgfrq[valcnt][0] = np.arccos(1-x)
            qgfrq[valcnt][1] = k
            indx = qgcnt - 1 - valcnt
            qgfrq[indx] = 2*np.pi - qgfrq[valcnt][0]            
            qgfrq[indx][1] = k
            valcnt += 1
       
# Get the number of edges
    Nedge = np.shape(Edges)[0] 
   
# Make an array for the qg eigenvectors coming from combinatorial ones.    
# For each vector there are Nedge edges.
# The edges are oriented from IntEdge[.][0] to IntEdge[.][1]
# On each edge the qg eigenvectors are a linear combination of 
# cos(omega x) and sin(omega x) or exp(iomega x) and exp(-iomega x)
# The coefficients are computed from the end vertex values a,b
# The first value is the cosine coefficient.  The second is the sine coefficient. 
    
  
    QGfrq = np.zeros(qgcnt)
    QGevec0 = 1j*np.zeros((Nedge,qgcnt,2))
    for j in np.arange(qgcnt):
        omega = qgfrq[j][0]
        QGfrq[j] = omega
        evecpnt = int(qgfrq[j][1])  
        for k in np.arange(Nedge):
            vert1 = int(Edges[k][0]) 
            vert2 = int(Edges[k][1])
            a = eigenvectors[vert1][evecpnt].copy()
            QGevec0[k][j][0] = a 
            b = eigenvectors[vert2][evecpnt].copy()
            QGevec0[k][j][1] =  (b - a*np.cos(omega))/(np.sin(omega))
           
 # Change from trigonometric coefficients to exponential coefficients    
    QGevec = 1j*np.zeros((Nedge,qgcnt,2))
    for j in np.arange(qgcnt):
       for k in np.arange(Nedge):
           QGevec[k][j][0] = 0.5*(QGevec0[k][j][0] - 1j*QGevec0[k][j][1])
           QGevec[k][j][1] = 0.5*(QGevec0[k][j][0] + 1j*QGevec0[k][j][1])    
                
    return(QGfrq,QGevec)

