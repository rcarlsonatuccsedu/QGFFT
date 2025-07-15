#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 16:04:40 2024

@author: robertcarlson
"""

import numpy as np
from scipy.linalg import null_space 

# Find the eigenfunctions with eigenvalues 0, pi^2 or (2pi)^2

def PiEqns(freq,Adj,Edges):
    
    Nedge = np.shape(Edges)[0]
        
# Find the degrees of the vertices 
    Nvert = int(np.shape(Adj)[0])   
    Deg = np.zeros(Nvert)
    Degsum = 0
    Degmax = 2
    for j in np.arange(Nvert):
        for k in np.arange(Nvert):
            Deg[j] = Deg[j] + Adj[j][k]
        if Deg[j] > Degmax:
            Degmax = int(Deg[j])
        Degsum = Degsum + Deg[j]
 
    
    # For each vertex, record its incident edges

    vertinfo = Edgerecord(Nvert,Degmax,Edges)
    
 
    # Construct the matrix encoding the vertex conditions.
        
    VertCond = np.zeros((int(Degsum),int(2*Nedge)))
 #   print('Degsum = ',int(Degsum))
    eqnpntr = 0
    # For every vertex, encode the continuity and derivative equations
    # For each vertex j there are Deg[j]-1 continuity conditions
    # and one derivative condition.
    
    # First the continuity conditions

    for j in np.arange(Nvert):    
        e1 = findedgenum(Edges,j,vertinfo[j][0][0])
        for k in np.arange(1,int(Deg[j])):
            VertCond[eqnpntr][int(e1)] = -1
            e2 = findedgenum(Edges,j,vertinfo[j][k][0])  
            if vertinfo[j][k][1] == vertinfo[j][0][1]:
                VertCond[eqnpntr][int(e2)] = 1
            else:     
                VertCond[eqnpntr][int(e2)] = np.cos(freq)
            eqnpntr += 1    
   
    # Then the derivative conditions
    if freq == np.pi:
        for j in np.arange(Nvert):
            e1 = int(findedgenum(Edges,j,vertinfo[j][0][0]))
            VertCond[eqnpntr][e1 + Nedge ] = 1
            for k in np.arange(1,int(Deg[j])):
                e2 = int(findedgenum(Edges,j,vertinfo[j][k][0]))
                VertCond[eqnpntr][e2 + Nedge] = 1
            eqnpntr += 1
            
    if freq == 2*np.pi:
        for j in np.arange(Nvert):
            e1 = int(findedgenum(Edges,j,vertinfo[j][0][0]))
            if vertinfo[j][0][1] == 1:
                VertCond[eqnpntr][e1 + Nedge ] = 1
            else:
                VertCond[eqnpntr][e1 + Nedge ] = -1
            for k in np.arange(1,int(Deg[j])):
                e2 = int(findedgenum(Edges,j,vertinfo[j][k][0]))
                if vertinfo[j][k][1] == vertinfo[j][0][1]:
                    VertCond[eqnpntr][e2 + Nedge] = VertCond[eqnpntr][e1 + Nedge ]
                else:                         
                    VertCond[eqnpntr][e2 + Nedge] = -VertCond[eqnpntr][e1 + Nedge ]
            eqnpntr += 1

    Nullsp = null_space(VertCond)
    
    Nullnum = np.shape(Nullsp)[1]
    
    # The output for each of the two freq values will be an array
    # with first index the edge number, the second index the vector number
    # and the third index indicating the (cosine,sine) coefficient
  
   
    Eigvecs = np.zeros((Nedge,Nullnum,2))
    for m in np.arange(Nullnum):
        for n in np.arange(Nedge):
            Eigvecs[n][m][0] = Nullsp[n][m] 
            Eigvecs[n][m][1] = Nullsp[n+Nedge][m]
                
    return(Eigvecs)

def Edgerecord(Nvert,Degmax,Edges):
    # vertinfo will record the edges incident on the vertices
    # The second component will be an adjacent vertex index
    # The last component is an orientation indicator
    NEdge = np.shape(Edges)[0]
    
    vertinfo = np.zeros((Nvert,Degmax,2))
    for j in np.arange(Nvert):
        epnt = 0
        for k in np.arange(NEdge):
            if Edges[k][0] == j:
                vertinfo[j][epnt][0] = Edges[k][1]
                vertinfo[j][epnt][1] = 1
                epnt += 1
            if Edges[k][1] == j:
                vertinfo[j][epnt][0] = Edges[k][0]
                vertinfo[j][epnt][1] = 2
                epnt += 1
        
    return(vertinfo)

# Find the index of an edge from the vertex pair m,n
def findedgenum(Edges,m,n):
    NEdge = np.shape(Edges)[0]
    for k in np.arange(NEdge):
        m1 = Edges[k][0]
        n1 = Edges[k][1]
        if (m == m1) and (n == n1):
            j=k
        if (m == n1) and (n == m1):
            j = k
    return(j)