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

# Construct the coefficient function sigma
def Phi(le):
    
    Nsamp = Constants.NSAMP
    twopi = Constants.TWOPI
    Nplus1 = Nsamp +  1
    phie = np.zeros(Nplus1)
    for n in np.arange(Nplus1):
        z = le*n/Nsamp
        phie[n] = z/le -(1-le)*np.sin(twopi*z/le)/twopi
    phiinv = Interp(phie)  
        
    return(phie,phiinv)

 
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
         wt = (m - low)/(high-low)
         phiinv[n] = wt*(lowind+1) + (1-wt)*lowind
         if (wt <= .5) and (lowind > 0):
             phiinv[n] = (lowind-1)*wt*(wt-1)/2 - lowind*(wt+1)*(wt-1) + (lowind+1)*wt*(wt+1)/2 
         elif (wt >= .5) and (lowind < (Nsamp-1)):
             wt = (m - high)/(high-low)
             phiinv[n] = lowind*wt*(wt-1)/2 - (lowind+1)*(wt+1)*(wt-1) + (lowind+2)*wt*(wt+1)/2 
     
     return(phiinv)

    