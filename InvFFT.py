#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
""" 
import numpy as np

# The input FT has size (Neval,Nsamp/2)
# After zero padding FT, an Nsamp size FFT is computed
def InvFT(oddcase,Nsamp,FT,Frqs,Evecs):
    
    Nfilt = int(Nsamp/2)
    Nfft = Nsamp
    Nfrq = np.shape(Frqs)[0]
    Nedge = np.shape(Evecs)[0]
    Nplus1 = Nsamp + 1

    IFT = 1j*np.zeros((Nedge,Nplus1))   
    IFTD = 1j*np.zeros((Nedge,Nplus1)) 
    IFTm = 1j*np.zeros((Nedge,Nplus1))  
    IFT1 = 1j*np.zeros((Nedge,Nfft))   
     
    ModFunc = 1j*np.zeros((Nfrq,Nfft))        
    
# Rescale oddcase
    FT[int(oddcase)][Nfilt-1] = np.sqrt(0.5)*FT[int(oddcase)][Nfilt-1]       
    
    for m in np.arange(Nfrq):
        for j in np.arange(Nfilt):              
            ModFunc[m][j] = FT[m][j]   
        for n in np.arange(Nedge):     
            IFT1[n][:] = np.fft.fft(ModFunc[m][:])  
        omega = Frqs[m]
        for n in np.arange(Nedge):
            Edata1 = Evecs[n][m][0]
            Edata2 = Evecs[n][m][1]
            for k in np.arange(1,Nsamp):  
                x = 1j*omega*k/Nsamp 
                z1 = np.exp(x) 
                z2 = np.exp(-x)
                IFTm[n][k] = z1*Edata1*IFT1[n][Nsamp-k] + z2*Edata2*IFT1[n][k]
            IFTm[n][0] = Edata1*IFT1[n][0] + Edata2*IFT1[n][0]
            x = 1j*omega
            z1 = np.exp(x) 
            z2 = np.exp(-x)
            IFTm[n][Nsamp] = z1*Edata1*IFT1[n][0] + z2*Edata2*IFT1[n][0]
            
        IFT = IFT + IFTm

    return (IFT)



