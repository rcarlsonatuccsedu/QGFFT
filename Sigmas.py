#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
import Plots, Constants
import matplotlib.pyplot as plt

# Construct the coefficient function sigma and related functions used for the
# transformed Laplacian coefficients.
def Sigdef(Edges):
    
    
    twopi = Constants.TWOPI
    Nsamp = Constants.NSAMP
    Nplus1 = Nsamp + 1
    Nedge = np.shape(Edges)[0]
    
    sigma = np.zeros((Nedge,Nplus1))
    sigmapr = np.zeros((Nedge,Nplus1))
    sigprpr = np.zeros((Nedge,Nplus1))
    sigma2 = np.zeros((Nedge,Nplus1))
    sigma3 = np.zeros((Nedge,Nplus1))
       
    for j in np.arange(Nedge):
        phie = np.zeros(Nplus1)
        le = Edges[j][2]
        fac = (1-le)/le 
        for n in np.arange(Nplus1):
            x = le*n/Nsamp
            # Both phie and phie inverse are arrays of length Nsamp + 1 
            phie[n] = x/le - (1 - le)*np.sin(twopi*x/le)/twopi  
        phiinv = Interp(phie)  
        for n in np.arange(Nplus1):
            z = phiinv[n]*le/Nsamp  # value of z = phi ^{-1}(x) 
            sigma[j][n] = 1/le - fac*np.cos(twopi*z/le)
            phi2 = twopi*(1-le)*np.sin(twopi*z/le)/(le**2) # phi''
            sigmapr[j][n] = phi2/sigma[j][n]  # sigma'
            phi3 = ((twopi**2)*(1-le)/(le**3))*np.cos(twopi*z/le) # phi '''
            sigprpr[j][n] = (phi3 - phi2*sigmapr[j][n])/(sigma[j][n]**2) 
    
    sigrt = np.sqrt(sigma)  # sqrt(sigma)
    sigsq = sigma*sigma  # sigma squared
    sigpot = sigmapr*sigmapr/4 - sigma*sigprpr/2  # potential term 

    return (sigrt,sigsq,sigpot)

# The function Interp will construct the inverse of phie
def Interp(phie):
    
    Nsamp = Constants.NSAMP
    Nplus1 = Nsamp +  1
    phiinv = np.zeros(Nplus1)
    
    low = phie[0]
    high = phie[1]
    lowind = int(0)
    for n in np.arange(Nplus1):
        m = n/Nsamp
        while m > high:
            lowind = lowind + 1
            low = phie[lowind]
            high = phie[lowind+1]        
        wt1 = (m - low)/(high-low)
        phiinv[n] = wt1*(lowind+1) + (1-wt1)*lowind
    
    
        
    return(phiinv)