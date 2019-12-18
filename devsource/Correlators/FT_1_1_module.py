#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 12:10:36 2019

@author: ppxsw1
"""

import numpy as np

def ft_1_1(phi, dx):
    Nx = phi.shape[0]
    Np = Nx
    i = complex(0,1)
    phi_p = np.zeros(phi.shape)
    for q in np.arange(Np):
        p = (2*(1 - np.cos(q*2*np.pi/(dx*Nx))))**2
        p /= dx**2
        p = np.sqrt(p)
        for j in np.arange(Nx):
            phi_p[q, :] += phi[j, :]*np.exp(i*j*dx*p)
    phi_p *= dx
    return phi_p

def ft(phi, p, dx):
    phi_p = np.zeros(phi.shape[1], dtype = complex)
    i = complex(0, 1)
    for j in np.arange(phi.shape[0]):
        phi_p += phi[j, :]*np.exp(i*j*dx*p)
    return phi_p

def calc_ddF(F, dt):
    dF = np.gradient(F, dt, axis = 0)
    ddF = np.gradient(dF, dt, axis = 1)
    ret = np.diag(ddF)
    return ret
