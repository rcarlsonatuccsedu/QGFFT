#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 08:42:52 2023

@author: robertcarlson
"""
import numpy as np


def MakeEdgeList(Adjac):

# Input is the initial adjacency matrix with nonnegative integer entries
# The output is a list of edges.
# Each edge is encoded as a pair of vertices, with smallest vertex index first.    
# In addition there is a list of original edge lengths.

# Count edges in the input graph description
# NAdjac is the number of original graph vertices   
# NEdge0 is number of edges in the original graph 
    NAdjac = np.shape(Adjac)[0]
    NEdge0 = 0
    for m in np.arange(NAdjac):
        for n in np.arange(m,NAdjac):
            if Adjac[m][n] > 0: NEdge0 += 1  
            
# Edge0 is an array of directed original edges
# Data is (vertex1, vertex2,length)
# Each edge is recorded as a pair of vertices, smaller vertex number first
# Elength is an array of the edge lengths
    Edge0 = np.zeros((NEdge0,3))    
    Elength = np.zeros(NEdge0)
    eIndx = 0
    for m in np.arange(NAdjac):
        for n in np.arange(m,NAdjac):
            if Adjac[m][n] > 0:
                Edge0[eIndx][0] = m
                Edge0[eIndx][1] = n
                Edge0[eIndx][2] = Adjac[m][n]
                eIndx += 1
    print('Number of edges = ',NEdge0)


    print('Edge number,vert1, vert2, length ')
    for kk in np.arange(NEdge0):
        k = int(kk)        
        print('Edge: ',k,Edge0[k][0],Edge0[k][1],Edge0[k][2])

    return(Edge0)

