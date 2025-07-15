#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
import Plots, InvFFT,FFT
import Constants, ChofV
import matplotlib.pyplot as plt

# Generate data on graph
def SchrdData(sigrt,Edges,Nedge,Frqs,Evecs):
    
    twopi = Constants.TWOPI
    Nsamp = Constants.NSAMP
    Nplus1 = int(Nsamp + 1)
    Nfft = int(Nsamp/2)
    
#Finit is the initial data     
    Finit = 1j*np.zeros((Nedge,Nplus1))
   
    c = Constants.PACKFRQ
    scale = .4
    frq = Constants.FRQ
    for k in np.arange(Nplus1):
        x = k/Nsamp
        Finit[3][k] = .4*(1-np.cos(twopi*x))*np.exp(c*x)
           
    return (Finit)

