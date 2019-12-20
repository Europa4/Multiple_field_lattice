# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 19:36:03 2019

@author: Woody
"""

import numpy as np
import pandas as pd
from multiprocessing import Pool
import time

import prot
import FT_1_1_module as FT
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

def observable(phi_data, header, p):
    #This is the observable you want to jam into the equation 26 from the below arxiv link
    i_phi = np.reshape(phi_data, (prot.Nx, prot.Nrp))
    i_phi = i_phi[:, :(prot.Nt - 2)]
    phi_k = FT.ft(i_phi, p, prot.deltaX)
    phi_nk = FT.ft(i_phi, -1*p, prot.deltaX)
    phi_k = np.reshape(phi_k, (1, prot.Nt - 2))
    phi_nk = np.reshape(phi_nk, (prot.Nt - 2, 1))
    propogator = np.matmul(phi_nk, phi_k)
    return propogator

def calc_expectation_and_error(n):
    expectation_observable = np.zeros((prot.Nt - 2, prot.Nt - 2), dtype = complex)
    file_error = np.copy(expectation_observable)
    #data location
     
    file_name = location + 'phi_' + str(n)
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

    number_of_iterations = phi_data.shape[0]
    #equation 26 in arxiv 1704.06404, |det J(0)| is ingored as it'll cancel between the numerator and the denominator
    phi_tilde = np.exp(aux_data[:, 2] + prot.j*aux_data[:, 3] - prot.j*aux_data[:, 1])

    denominator = np.mean(phi_tilde)
    numerator = np.zeros((number_of_iterations, prot.Nt-2, prot.Nt-2), dtype = complex)
    p = 0
    for i in np.arange(int(number_of_iterations)):
        Obs = observable(phi_data[i, :], header, p)
        numerator[i, :] = Obs*phi_tilde[i] #this is the unaveraged numerator
    #This calculates the expression fully for one initiatlisation
    expectation_observable = np.mean(numerator, axis = 0)/denominator

    #error analysis for this file
    for i in np.arange(prot.Nt - 2):
        for j in np.arange(prot.Nt - 2):
           file_error[i, j] = jackknife(np.real(numerator[:, i, j]/denominator), jackknife_block_length) + prot.j*jackknife(np.imag(numerator[:, i, j]/denominator), jackknife_block_length)
    return [expectation_observable, file_error]

n_files = 100
jackknife_block_length = 250
expectation_observable = np.zeros((n_files, prot.Nt - 2, prot.Nt - 2), dtype = complex)
file_error = np.zeros((n_files, prot.Nt - 2, prot.Nt - 2), dtype = complex)
x_range = np.arange(prot.Nt) + 1
n_threads = 2
p = Pool(n_threads)
output = p.map(calc_expectation_and_error, np.arange(n_files))

for i in np.arange(n_files):
    expectation_observable[i, :] = output[i][0]
    file_error[i, :] = output[i][1]

final_expectation = np.mean(expectation_observable, axis = 0)

final_error = np.sqrt(np.sum(np.real(file_error)**2, axis = 0) + prot.j*np.sum(np.imag(file_error)**2, axis = 0))/np.sqrt(n_files)

F = np.diag(final_expectation)
ddF = FT.calc_ddF(final_expectation, prot.deltaT)

w_squared = ddF/F
c = np.sqrt(1 - 0.25*(prot.deltaT)**2*w_squared)
n = c*np.sqrt(ddF*F) - 0.5

print(n)
