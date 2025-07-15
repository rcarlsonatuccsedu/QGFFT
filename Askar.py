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
def E1(Frqs,Evecs,k,thyme,sigsq,sigpot,Soln):
    
    Nsamp = Constants.NSAMP
    DiffCoef = Constants.COEF
   
    # First term
    D1 = NewLap.Derive2(Frqs,Evecs,Soln[k-1][:][:])
    T1 = sigsq*D1
    
# Second term
  
    T2 = Soln[k-1][:][:]*sigpot
    
    V1 = DiffCoef*(T1 + T2)
  
    if k == 1:
        D11 = NewLap.Derive2(Frqs,Evecs,Soln[k-1][:][:]+thyme*V1)
        T11 = sigsq*D11
        T22 = sigpot*(Soln[k-1][:][:]+thyme*V1)
        V2 = DiffCoef*(T11 + T22)
        Soln[k][:][:] = Soln[k-1][:][:] + thyme*(V1+V2)/2
    else:
        Soln[k][:][:] = Soln[k-2][:][:] + 2*thyme*V1
    
    
    return(Soln)