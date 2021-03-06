# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:35:31 2020

@author: Woody
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 15:46:50 2020

@author: ppxsw1
"""

import numpy as np
import pandas as pd
from multiprocessing import Pool
import matplotlib.pyplot as plt

import prot
#location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/lambda_0/Tau_1-50/Delta_0-1500/'
#location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/Data_1_1/'
location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/long/'
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
    #currently produces a matrix of <phi+_i phi+_j> correlators
    phi = phi_data[:, 0:(prot.Nt - 2)] #taking only the upper branch of each layer
    field_0 = np.reshape(header[1:(prot.Nx + 1)], (prot.Nx, 1))
    field_1 = np.reshape(header[(2*prot.Nx + 1):(3*prot.Nx + 1)], (prot.Nx, 1))
    phi = np.concatenate((field_0, field_1, phi), axis = 1)
    phi = np.reshape(phi, (phi.size, 1))
    phi_t = np.transpose(phi)
    correlator = np.matmul(phi, phi_t)
    return correlator

def calc_expectation_and_error(n_file):
     
    file_name = location + 'phi_' + str(n_file)
    #This is all just using pandas to load the CSV files, it's completely compatable with the outputs from the C++ code
    #It outputs it as a triplet of numpy arrays, header, aux_data, and phi_data
    data = pd.read_csv(file_name, header = None, skiprows = 1)
    header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
    data = np.array(data)
    aux_data = data[:, 2*prot.Nrt:]
    phi_data = data[:, :2*prot.Nrt]
    #This gets it from the CSV state of "real, imag, real, imag" to "real + i*imag, real + i*imag"
    temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
    phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]

    #need to find the number of MC rows in the file
    number_of_iterations = phi_data.shape[0]
    #equation 26 in arxiv 1704.06404, |det J(0)| is ingored as it'll cancel between the numerator and the denominator
    phi_tilde = np.exp(aux_data[:, 2] + prot.j*aux_data[:, 3] - prot.j*aux_data[:, 1])
    denominator = np.mean(phi_tilde)
    numerator = np.zeros((number_of_iterations, prot.Nx*prot.Nt, prot.Nx*prot.Nt), dtype = complex)
    
    for i in np.arange(number_of_iterations):
        step_phi_data = np.reshape(phi_data[i, :], (prot.Nx, prot.Nrp))
        F = observable(step_phi_data, header)
        numerator[i, :, :] = F*phi_tilde[i]
        
    averaged_n = np.mean(numerator, axis = 0)/denominator
    print("file", str(n_file), "analysed")
    return averaged_n

n_files = 5
jackknife_block_length = 250
expectation_observable = np.zeros((n_files, prot.Nx*prot.Nt, prot.Nx*prot.Nt), dtype = complex)
file_error = np.zeros((n_files, prot.Nx, prot.Nt), dtype = complex)
x_range = (np.arange(prot.Nt) + 1)*prot.deltaT
n_threads = 1
p = Pool(n_threads)
output = p.map(calc_expectation_and_error, np.arange(n_files))
#for i in np.arange(n_files):
#    calc_expectation_and_error(i)
#    print("finished iteration", str(i))

for i in np.arange(n_files):
    expectation_observable[i, :, :] = output[i]
final_expectation = np.mean(expectation_observable, axis = 0)
np.save("spacial_correlator", expectation_observable)

plt.figure(0)
plt.plot(np.real(final_expectation[0, :prot.Nt]))
plt.plot(np.imag(final_expectation[0, :prot.Nt]))
plt.savefig("test.png")
