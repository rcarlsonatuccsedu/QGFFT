#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 08:42:52 2023

@author: robertcarlson
"""
import numpy as np
# Input: list of edges given as vertex pairs
# Output: Adjacency matrix, list of combinatorial Laplacian eigenvalues, 
#list of combinatorial Laplacian eigenvectors
def CombEval(Edges):
# Compute the eigenvalues of the combinatorial Laplacian    
# starting from the edge list of the integerized graph    
# The number of vertices is one more than the maximum vertex number. 

    NEdge = np.shape(Edges)[0]
    NVert = int(np.max(Edges) + 1)

# Make the adjacency matrix
    Adj = np.zeros((NVert,NVert))
    for k in np.arange(NEdge):
        m = int(Edges[k][0])
        n = int(Edges[k][1])
        Adj[m][n] = 1
        Adj[n][m] = 1
  
# Denote the degree operator by T    
# Compute the vector of vertex degrees and the diagonal of T^{-1/2}  
    T0 = np.zeros(NVert)
    T1 = np.zeros(NVert)
    for j in np.arange(NVert):
        for k in np.arange(NVert): 
            T1[j] = T1[j] + Adj[j][k]
        if T1[j] > 0:
            T0[j] = 1.0/np.sqrt(T1[j])

# Make the combinatorial Laplacian I - T^{-1/2}AT^{-1/2}    
# which is a real symmetric matrix
    L1 = np.zeros((NVert,NVert))
    for j in np.arange(NVert):  
        for k in np.arange(NVert):
           L1[j][k] = T0[j]*Adj[j][k]*T0[k]
    Id = np.zeros((NVert,NVert))
    for j in np.arange(NVert):
      Id[j][j] = 1
    Lap = Id - L1
    
# Get eigenvalues and (dot product) orthonormal eigenvectors of Lap
    eigenvalues,eigenvectors = np.linalg.eigh(Lap)

# Get eigenvectors for matrix I - T^{-1}A
# These eigenvectors are orthonormal for the degree weighted inner product
    Combevecs =   1j*np.zeros((NVert,NVert))   
    for j in np.arange(NVert):
        for k in np.arange(NVert):
            Combevecs[j][k] = T0[j]*eigenvectors[j][k]

    return(Adj,eigenvalues,Combevecs)


