#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 08:42:52 2019

@author: robertcarlson
"""
import numpy as np
NSAMP = 32
TWOPI = 2*np.pi
COEF = 1j
# The 1/16 matches the fraction of a cycle for eigenfunction with 
# eigenvalue (pi)^2 
#TIMESCALE = 1/(16*np.pi)
TIMESCALE = .006
FRQ = 1
LE = np.sqrt(2)/2
PACKFRQ = 20j*TWOPI*(NSAMP/32)
LAM = COEF*((FRQ*TWOPI/LE)**2)
LAM1 = COEF*((FRQ*np.pi)**2)