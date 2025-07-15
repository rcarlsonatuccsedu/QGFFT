#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
""" 
import numpy as np
import matplotlib.pyplot as plt

def graphFT(Nsamp,Func,PiEvals,PiEvecs):
    
    Nfft = Nsamp
    Neval = np.shape(PiEvals)[0]
    Nedge = np.shape(PiEvecs)[0]

    print('Nedge = ',Nedge)

    FT = 1j*np.zeros((Neval,Nfft))   
    FT1 = 1j*np.zeros((Neval,Nfft)) 
    FT2 = 1j*np.zeros((Neval,Nfft))
    
  # Filter 0 has no higher frequency modes
    Filt0 = 0
    for j in np.arange(Nedge):        
      for k in np.arange(Nsamp+1):
          wt = 1
          if (k == 0) or (k == Nsamp):
               wt = .5   
          Filt0 = Filt0 + wt*Func[j][k]
            
    for m in np.arange(1,Neval):
        ModFunc1 = 1j*np.zeros((Nedge,Nsamp))        
        ModFunc2 = 1j*np.zeros((Nedge,Nsamp))
        omega = PiEvals[m]
        if (m == 11):
            print('omega =',omega)
        for n in np.arange(Nedge):
            Edata1 = np.conj(PiEvecs[n][m][0])
            Edata2 = np.conj(PiEvecs[n][m][1])
            if (m == 11):
                print('Edata = ',Edata1,Edata2)
            z1 = np.cos(omega) - 1j*np.sin(omega)
            z2 = np.cos(omega) + 1j*np.sin(omega)
            ModFunc1[n][0] = .5*Edata1*(Func[n][0] + z1*Func[n][Nsamp])  
            ModFunc2[n][0] = .5*Edata2*(Func[n][0] + z2*Func[n][Nsamp])  
            for k in np.arange(1,Nsamp):      
               x = omega*k/Nfft
               z1 = np.cos(x) - 1j*np.sin(x)
               z2 = np.cos(x) + 1j*np.sin(x)
               ModFunc1[n][k] = Edata1*z1*Func[n][k]   
               ModFunc2[n][k] = Edata2*z2*Func[n][k]   
            FT1[m][:] = FT1[m][:] + np.fft.fft(ModFunc1[n][:])  
            FT2[m][:] = FT2[m][:] + np.fft.fft(ModFunc2[n][:])   
    
    FT[0][0] = Filt0
    for m in np.arange(1,Neval):
        FT[m][0] = FT1[m][0] + FT2[m][0]
        for k in np.arange(1,int(Nfft/2)):
            FT[m][k] = FT1[m][k] + FT2[m][Nfft - k]   
            

    return (FT)



def DisplayFT(eigval,PiFT):

    np.set_printoptions(3)
    
    Neval = np.shape(PiFT)[0]
    Nsamp = np.shape(PiFT)[1]
    
    xcoord = np.arange(Nsamp)     
    plt.figure()
    plt.plot(xcoord,np.real(PiFT[eigval][:]),np.imag(PiFT[eigval][:]))

    return ()