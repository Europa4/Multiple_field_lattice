#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  5 16:53:03 2018

@author: ppxsw1
"""
import numpy as np
Nx = 1
Nt = 10
Nb = 2
Npath = 2*Nt
Ntot = Nx*Npath
Nrp = Npath - 4
Nrt = Nx*Nrp
deltaT = 0.75
deltaX = 1.0
Lambda = 0.0
m = 1.0
delta = 0.1
j = complex(0,1)
dt = deltaT*np.concatenate((np.ones(Nt-2), -np.ones(Nt-2)))
