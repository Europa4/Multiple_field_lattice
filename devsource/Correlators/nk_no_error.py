# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 08:45:06 2020

@author: Woody
"""

import numpy as np
import matplotlib.pyplot as plt

import prot

unaveraged_data = np.load("expectation_observable.npy")
averaged_data = np.mean(unaveraged_data, axis = 0)

t_range = np.arange(prot.Nt)*prot.deltaT 

for i in np.arange(prot.Nx):
    plt.figure(i)
    plt.plot(t_range, np.real(averaged_data[i, :]), 'x', label = 'real')
    plt.plot(t_range, np.imag(averaged_data[i, :]), 'x', label = 'imag')
    plt.title(r"$n_k$ for $k = $" + str(i))
    plt.xlabel("t")
    plt.ylabel(r"$n_k$")
    plt.legend()
    plt.savefig("no_error_n_" + str(i))