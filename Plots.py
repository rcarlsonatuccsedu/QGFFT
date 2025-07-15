#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
import matplotlib.pyplot as plt
import Constants

def Loopy(Eorder,Edir,Esamp,NewSoln):
    
    Tsteps = np.shape(NewSoln)[0]
    Nedge = np.shape(Esamp)[0]
    Nsamp = Constants.NSAMP

    print('Loopy Tsteps:',Tsteps)

    Npnt = 0
    for j in np.arange(Nedge):
        Npnt = Npnt + Esamp[j] + 1
    Npnt = int(Npnt)    

    buf = 1j*np.zeros(Npnt) 
 
    w1 = np.zeros(Npnt)
    Np1 = Nsamp + 1 
    for k in np.arange(4*Np1,6*Np1):
        w1[k] = -.05
    pt1 = 9*Np1 + Esamp[5] + Esamp[8]
    for k in np.arange(9*Np1,pt1):
        w1[k] = -.05
    
    for j in np.arange(0,Tsteps,20):
        pntr = 0        
        for k in np.arange(Nedge):
            kk = Eorder[k]
            N = Esamp[kk]       
            for m in np.arange(N+1):                
                buf[pntr] = NewSoln[j][kk][m]
                if int(Edir[k]) == int(1):
                    mk = int(Esamp[kk]) - m 
                    buf[pntr] = NewSoln[j][kk][mk]
                pntr = pntr + 1
    
        
        x = np.arange(Npnt)
        w = np.zeros(Npnt)
        for j1 in np.arange(Npnt):
            w[j1] = np.abs(buf[j1])**2     
        
       
        plt.figure()
        plt.ylim(-.1,1)
        plt.title("Packet hits asymmetric box")   
        plt.plot(x,w,'k',x,w1,'k')
        plt.show()
        print('Tsteps = ',j)
        input('Hit return')
        
    return()


