#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 14:47:53 2019

@author: ppxsw1
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import prot

def jackknife(data, block_length):
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
    #This is the observable you want to jam into the equation 26 from the above arxiv link
    q_valid_range = prot.Nt - 2
    quantum_quantum_vector = phi_data[1:q_valid_range] - np.flip(phi_data[(q_valid_range + 1):])
    quantum_quantum_vector = np.concatenate(([0, 2*phi_data[0]],quantum_quantum_vector, [0]))
    
    classical_vector = 0.5*(phi_data[1:q_valid_range] + np.flip(phi_data[(q_valid_range + 1):]))
    classical_vector = np.concatenate(([header[1], header[3]], classical_vector, [phi_data[q_valid_range]]))
    
    classical_vector = np.reshape(classical_vector, [classical_vector.size, 1])
    quantum_quantum_vector = np.reshape(quantum_quantum_vector, [1, quantum_quantum_vector.size])
    classical_quantum_correlator = np.matmul(classical_vector, quantum_quantum_vector)
    
    quantum_classical_correlator = np.matmul(np.transpose(quantum_quantum_vector), np.transpose(classical_vector))
    
    classical_classical_correlator = np.matmul(classical_vector, np.transpose(classical_vector))
    forward_correlator = classical_classical_correlator + 0.5*(classical_quantum_correlator + quantum_classical_correlator)
    return forward_correlator


n_files = 10
jackknife_block_length = 250
expectation_observable = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
file_error = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
x_range = np.arange(prot.Nt) + 1
location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/lambda_0/Tau_1-50/Delta_0-1500/'
for n in np.arange(n_files):
    #data location
    file_name = location + 'phi_' + str(n)
    #This is all just using pandas to load the CSV files, it's completely compatable with the outputs from the C++ code
    #It outputs it as a triplet of numpy arrays, header, aux_data and phi_data
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
    numerator = np.zeros((number_of_iterations, prot.Nt, prot.Nt), dtype = complex)
    
    for i in np.arange(int(number_of_iterations)):
        Obs = observable(phi_data[i, :], header)
        numerator[i, :, :] = Obs*phi_tilde[i] #this is the unaveraged numerator
    #This calculates the expression fully for one initiatlisation
    expectation_observable[n, :, :] = np.mean(numerator, axis = 0)/denominator
    
    #error analysis for this file
    for i in np.arange(prot.Nt):
        for j in np.arange(prot.Nt):
            file_error[n, i, j] = jackknife(np.real(numerator[:, i, j]/denominator), jackknife_block_length) + prot.j*jackknife(np.imag(numerator[:, i, j]/denominator), jackknife_block_length)

#this now combines all the initialisations
final_expectation = np.mean(expectation_observable, axis = 0)

final_error = np.sqrt(np.sum(np.real(file_error)**2, axis = 0) + prot.j*np.sum(np.imag(file_error)**2, axis = 0))/np.sqrt(n_files)

#plotting, whatever is suitable for your observable
for i in np.arange(prot.Nt - 2) + 1:
    plt.figure(i)
    plt.errorbar(x_range, np.real(final_expectation[:, i]).flatten(), yerr = np.real(final_error[:, i]), label = "real", fmt = '', linestyle = 'None', marker = 'o', capsize = 5)
    plt.errorbar(x_range, np.imag(final_expectation[:, i]).flatten(), yerr = np.imag(final_error[:, i]), label = "imag", fmt = '', linestyle = 'None', marker = 'o', capsize = 5)
    plt.axhline(0,color='black',linewidth = 0.5)
    plt.title('classical quantum i = ' + str(x_range[i]))
    plt.xlabel('j')
    plt.ylabel(r'$\langle \phi^{cl}_j \phi^{q}_i \rangle$')
    plt.legend()    
    plt.savefig(location + 'bf_forward_' + str(i))
    plt.plot()


