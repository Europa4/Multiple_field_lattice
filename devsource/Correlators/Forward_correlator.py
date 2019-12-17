#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 16:50:35 2019

@author: ppxsw1
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:47:53 2019

@author: ppxsw1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from multiprocessing import Pool

import prot
location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/lambda_0/Tau_1-50/Delta_0-1500/'

def jackknife(data, block_length = 250):
    #takes a 1D array of real data (for example all the [i,j]th correlators from a single initialisation) and performs jackknife analysis on the mean
    number_of_blocks = int(data.size/block_length)
    mu = np.mean(data)
    blocked_data = np.reshape(data, (number_of_blocks, block_length))
    reduced_mu = np.zeros(number_of_blocks)
    for i in np.arange(number_of_blocks):
        reduced_mu[i] = np.mean(np.delete(blocked_data, i, axis = 0))
    error = np.sqrt((number_of_blocks - 1)*np.sum((reduced_mu - mu)**2)/number_of_blocks)
    return error

def observable(phi_data, header):
    #This is the observable you want to jam into the equation 26 from the below arxiv link
    #currently set up to calculate the feynman propogator from the zeroth site
    phi_0 = header[1]
    phi_1 = header[3]
    #propogator = np.concatenate(([phi_0**2, phi_0*phi_1], phi_0*phi_data[1:int(prot.Nrp/2 + 1)]))
    propogator = np.concatenate(([phi_0*phi_data[1], phi_1*phi_data[1]], phi_data[1]*phi_data[1:int(prot.Nrp/2 + 1)]))
    return propogator

def calc_expectation_and_error(n):
    expectation_observable = np.zeros(prot.Nt, dtype = complex)
    file_error = np.copy(expectation_observable)
    #data location
    file_name = location + 'phi_' + str(n)
    #This is all just using pandas to load the CSV files, it's completely compatable with the outputs from the C++ code
    #It outputs it as a triplet of numpy arrays, header, aux_data, and phi_data
    data = pd.read_csv(file_name, header = None, skiprows = 1)
    header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
    data = np.array(data)
    aux_data = data[:, 2*prot.Nrp:]
    phi_data = data[:, :2*prot.Nrp]
    #This gets it from the CSV state of "real, imag, real, imag" to "real + i*imag, real + i*imag"
    temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
    phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]
    number_of_iterations = phi_data.shape[0]
    
    #equation 26 in arxiv 1704.06404, |det J(0)| is ingored as it'll cancel between the numerator and the denominator
    phi_tilde = np.exp(aux_data[:, 2] + prot.j*aux_data[:, 3] - prot.j*aux_data[:, 1])
    denominator = np.mean(phi_tilde)
    numerator = np.zeros((number_of_iterations, prot.Nt), dtype = complex)
    
    for i in np.arange(int(number_of_iterations)):
        Obs = observable(phi_data[i, :], header)
        numerator[i, :] = Obs*phi_tilde[i] #this is the unaveraged numerator
    #This calculates the expression fully for one initiatlisation
    expectation_observable = np.mean(numerator, axis = 0)/denominator
    
    #error analysis for this file
    for i in np.arange(prot.Nt):
        file_error[i] = jackknife(np.real(numerator[:, i]/denominator), jackknife_block_length) + prot.j*jackknife(np.imag(numerator[:, i]/denominator), jackknife_block_length)
    return [expectation_observable, file_error]

n_files = 10
jackknife_block_length = 250
expectation_observable = np.zeros((n_files, prot.Nt), dtype = complex)
file_error = np.zeros((n_files, prot.Nt), dtype = complex)
x_range = np.arange(prot.Nt) + 1
n_threads = 2
p = Pool(n_threads)
output = p.map(calc_expectation_and_error, np.arange(n_files))

for i in np.arange(n_files):
    expectation_observable[i, :] = output[i][0]
    file_error[i, :] = output[i][1]

final_expectation = np.mean(expectation_observable, axis = 0)

final_error = np.sqrt(np.sum(np.real(file_error)**2, axis = 0) + prot.j*np.sum(np.imag(file_error)**2, axis = 0))/np.sqrt(n_files)

#plotting, whatever is suitable for your observable
#np.savez('burges_dat.npz', final_expectation, final_error)
plt.figure(0)
plt.errorbar(0.75*x_range, np.real(final_expectation).flatten(), yerr = np.real(final_error).flatten(), label = "real", fmt = '', linestyle = 'None', marker = 'o', capsize = 5)
plt.errorbar(0.75*x_range, np.imag(final_expectation).flatten(), yerr = np.imag(final_error).flatten(), label = "imag", fmt = '', linestyle = 'None', marker = 'o', capsize = 5)
plt.axhline(0, color='black',linewidth = 0.5)
plt.title(r'Forward correlator $\langle \phi^0 \phi^t \rangle$')
plt.xlabel('t')
plt.ylabel(r'$\langle \phi_0 \phi_j \rangle$')
plt.legend()
plt.plot()
plt.savefig(location + 'forward_propogator')