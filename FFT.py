 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
""" 
import numpy as np
import matplotlib.pyplot as plt

# For each fundamental eigenvalue and each edge 
# two FFTs of size Nsamp are computed and summed over edges
# Half of the FFT outputs are combined into one
# QGFFT of size Nsamp/2
def QGFT(Nsamp,Finit,Frqs,Evecs):
    
    Nfft = int(Nsamp/2)
    Nfrq = np.shape(Frqs)[0]
    Nedge = np.shape(Evecs)[0]
    Nplus1 = Nsamp + 1

    FT = 1j*np.zeros((Nfrq,Nfft))   
    FTinit = 1j*np.zeros((Nfrq,Nfft)) 
    FT1 = 1j*np.zeros((Nfrq,Nsamp)) 
    FT2 = 1j*np.zeros((Nfrq,Nsamp))
    
    ModFunc1 = 1j*np.zeros((Nfrq,Nedge,Nsamp))        
    ModFunc2 = 1j*np.zeros((Nfrq,Nedge,Nsamp))    
    
    for m in np.arange(Nfrq):
        omega = Frqs[m]
        for n in np.arange(Nedge):
            Edata1 = np.conj(Evecs[n][m][0])
            Edata2 = np.conj(Evecs[n][m][1])
            # handle endpoints
            x = 1j*omega
            z1 = np.exp(-x)
            z2 = np.exp(x)
            ModFunc1[m][n][0] = .5*Edata1*(Finit[n][0] + z1*Finit[n][Nsamp])   
            ModFunc2[m][n][0] = .5*Edata2*(Finit[n][0] + z2*Finit[n][Nsamp])  
            for k in np.arange(1,Nsamp):         
               x = 1j*omega*k/Nsamp
               z1 = np.exp(-x)
               z2 = np.exp(x)
               ModFunc1[m][n][k] = Edata1*z1*Finit[n][k]   
               ModFunc2[m][n][k] = Edata2*z2*Finit[n][k]   
            FT1[m][:] = FT1[m][:] + np.fft.fft(ModFunc1[m][n][:])  
            FT2[m][:] = FT2[m][:] + np.fft.fft(ModFunc2[m][n][:])   
    
    for m in np.arange(Nfrq):
        FT[m][0] = FT1[m][0] + FT2[m][0]
        for k in np.arange(1,Nfft):
            FT[m][k] = FT1[m][k] + FT2[m][Nsamp - k]
            
    # The eigenvalue 0 has no higher frequency terms        
    for k in np.arange(1,Nfft):
        FT[0][k] = 0
               
    # Adjust the oddball term
    for m in np.arange(Nfrq):
        if abs(Frqs[m] - 2*np.pi) < 10**(-10): 
            if (abs(Evecs[0][m][0]) > 10**(-10)) and abs(Evecs[0][m][0] - Evecs[0][m][1]) < 10**(-10):
                oddcase = m
                FT[m][Nfft-1] = np.sqrt(0.5)*FT[m][Nfft-1]
                     
        
    return (oddcase,FT)

    


