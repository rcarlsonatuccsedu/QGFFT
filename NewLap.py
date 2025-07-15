#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 16:18:14 2024

@author: robertcarlson
"""
import numpy as np  
import FFT, InvFFT
import Plots, Constants
def Derive(Frqs,Evecs,Crrnt):
    
    twopi = Constants.TWOPI
    Nsamp = Constants.NSAMP
    Np1 = int(Nsamp+1)
    # We need data from Evecs
    Nedges = np.shape(Evecs)[0]
    frqcnt = np.shape(Evecs)[1]
    highfrqcnt = int(Nsamp/2)
    
    # In order to compute derivatives we'll use the QGFFT coefficient
    # of Finit
    # The Fourier coefficients appear in an array of size (Nfrq,Nsamp/2)
    oddcase,FT = FFT.QGFT(Nsamp,Crrnt,Frqs,Evecs)
    
    mshape = np.shape(FT) 
   
    # Compute coefficients times frequencies
    #Then take FFT
    FTA = 1j*np.zeros((mshape[0],Nsamp)) 
    FT1 = 1j*np.zeros((frqcnt,Nsamp))
    FT2 = 1j*np.zeros((frqcnt,Nsamp))
    FS1 = 1j*np.zeros((frqcnt,Np1))
    FS2 = 1j*np.zeros((frqcnt,Np1))
    
    for funfrqndx in np.arange(1,frqcnt):
        omega0 = Frqs[funfrqndx]
        for m in np.arange(highfrqcnt):
            omegakm = omega0 + twopi*m
            FTA[funfrqndx][m] = FT[funfrqndx][m]*omegakm  
        FT1[funfrqndx][:] = np.sqrt(Nsamp)*np.fft.ifft(FTA[funfrqndx][:])  
        FT2[funfrqndx][:] = (1/np.sqrt(Nsamp))*np.fft.fft(FTA[funfrqndx][:]) 
        for n in np.arange(Nsamp):
            FS1[funfrqndx][n] = FT1[funfrqndx][n]
        FS1[funfrqndx][Nsamp] = FT1[funfrqndx][0]    
        for n in np.arange(Nsamp):
            FS2[funfrqndx][n] = FT2[funfrqndx][n]
        FS2[funfrqndx][Nsamp] = FT2[funfrqndx][0]  
    
    Fprime = 1j*np.zeros((Nedges,Np1))
    for edgndx in np.arange(Nedges):    
        for funfrqndx in np.arange(1,frqcnt):
            g = 1j*Evecs[edgndx][funfrqndx][0]
            d = 1j*Evecs[edgndx][funfrqndx][1]
            x = np.arange(Np1)/Nsamp
            c1 = g*np.exp(1j*Frqs[funfrqndx]*x)
            c2 = d*np.exp(-1j*Frqs[funfrqndx]*x)
            Fprime[edgndx][:] = Fprime[edgndx][:] + c1*FS1[funfrqndx][:] - c2*FS2[funfrqndx][:] 
            
    return (Fprime)

def Derive2(Frqs,Evecs,Crrnt):
    
    twopi = Constants.TWOPI
    Nsamp = Constants.NSAMP
    # We need data from Evecs
    frqcnt = np.shape(Evecs)[1]
    highfrqcnt = int(Nsamp/2)
    
    oddcase,FT = FFT.QGFT(Nsamp,Crrnt,Frqs,Evecs)
    
    FTA = 1j*np.zeros((frqcnt,highfrqcnt)) 
    
    for funfrqndx in np.arange(frqcnt):
        omega0 = Frqs[funfrqndx]
        for m in np.arange(highfrqcnt):
            omega2 = (omega0 + twopi*m)*(omega0 + twopi*m)
            FTA[funfrqndx][m] = FT[funfrqndx][m]*omega2 
    
    Fprpr = InvFFT.InvFT(oddcase,Nsamp,FTA,Frqs,Evecs)
    
    return(Fprpr)