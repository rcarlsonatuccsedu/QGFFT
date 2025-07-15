#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
import Constants
import matplotlib.pyplot as plt

# Rescale the solution 
def SMap(sigrt,Soln):
    
    Nsamp = Constants.NSAMP
    Np1 = Nsamp+1
    Tsteps = np.shape(Soln)[0]
    Nedge = np.shape(Soln)[1]
    Np1 = Nsamp + 1
    ScSoln = 1j*np.zeros((Tsteps,Nedge,Np1))
    for k in np.arange(Tsteps):
        ScSoln[k][:][:] = Soln[k][:][:]/sigrt

    return(ScSoln)