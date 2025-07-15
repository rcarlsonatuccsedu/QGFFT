#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: robertcarlson
"""

import numpy as np
import EdgeConvert

# An edge weighted input graph is selected.
# A second copy with lengths $1$ is made.
# The second copy is used to construct a list of length 1 edges,
# each edge encoded as a pair of vertices, with smallest vertex index first. 

def firststuff():

    # Make a weighted adjacency matrix with positive edge lengths

    L = int(7)
    Nvert = int(L+3)
    Adjac = np.zeros((Nvert,Nvert))
    for k in np.arange(1,L+2):
        Adjac[k-1][k] = 1
    Adjac[3][9] = np.pi/4    
    Adjac[0][L+1] = 1       
    Adjac[5][L+2] = np.sqrt(2)/2 

    for j in np.arange(Nvert):
        for k in np.arange(Nvert): 
            if Adjac[j][k] > 0:
                Adjac[k][j] = Adjac[j][k]
               
    # Make second copy with edge lengths 1
    Nvert = np.shape(Adjac)[0]
    Adjac0 = np.zeros((Nvert,Nvert))

    for m in np.arange(Nvert):
        for n in np.arange(Nvert):
            if Adjac[m][n] > 0:
                Adjac0[m][n] = 1

  # Form an array describing the graph edges.
  # for the integerized graph
    Edges = EdgeConvert.MakeEdgeList(Adjac)      

  # Define the order in which edge data will be displayed when solutions are computed.
    Enum = 11
    Eorder = [1,0,2,3,4,6,7,9,10,5,8]
    Edir = np.zeros(Enum)   
  # Some edge displays may require reversal of edge direction.  
    Edir[10] = 1

    return(Edges,Adjac0,Eorder,Edir)