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

def observable(phi_data, header, field_number):
    #This is the observable you want to jam into the equation 26 from the above arxiv link
    #currently set up to calculate the classical classical correlator for the 0D case
    site_0 = 2*field_number + 1
    site_1 = 2*field_number + 5
    classical_vector = np.reshape(np.concatenate((0.5*(phi_data[1:int(prot.Nrp/2)] + np.flip(phi_data[int(prot.Nrp/2 + 1):], axis = 0)), [phi_data[int(prot.Nrp/2)]])), (1, int(prot.Nrp/2)))
    classical_vector = np.append([header[site_0], header[site_1]], classical_vector)
    classical_vector = np.reshape(classical_vector, (1, classical_vector.size))
    classical_vector_transpose = np.reshape(classical_vector, (int(prot.Nrp/2) + 2,1))
    classical_classical_correlator = np.matmul(classical_vector_transpose, classical_vector)
    return classical_classical_correlator

n_files = 50
jackknife_block_length = 250
phi_expectation_observable = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
chi_expectation_observable = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
phi_file_error = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
chi_file_error = np.zeros((n_files, prot.Nt, prot.Nt), dtype = complex)
x_range = np.arange(prot.Nt) + 1
location = '../Data/2_test/'

for n in np.arange(n_files):
    #data location
    file_name = location + 'phi_' + str(n)
    #This is all just using pandas to load the CSV files, it's completely compatable with the outputs from the C++ code
    #It outputs it as a triplet of numpy arrays, header, aux_data and phi_data
    data = pd.read_csv(file_name, header = None, skiprows = 1)
    header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
    data = np.array(data)
    #aux_data = data[:, 2*prot.Nrp:]
    aux_data = data[:, -4:]
    phi_data = data[:, :2*prot.Nrp]
    chi_data = data[:, 2*prot.Nrp:4*prot.Nrp]
    #This gets it from the CSV state of "real, imag, real, imag" to "real + i*imag, real + i*imag"
    temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
    phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]
    number_of_iterations = phi_data.shape[0]

    temp = np.reshape(chi_data, (chi_data.shape[0], int(chi_data.shape[1]/2), 2))
    chi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]

    
    #equation 26 in arxiv 1704.06404, |det J(0)| is ingored as it'll cancel between the numerator and the denominator
    phi_tilde = np.exp(aux_data[:, 2] + prot.j*aux_data[:, 3] - prot.j*aux_data[:, 1])
    denominator = np.mean(phi_tilde)
    phi_numerator = np.zeros((number_of_iterations, prot.Nt, prot.Nt), dtype = complex)
    chi_numerator = np.zeros((number_of_iterations, prot.Nt, prot.Nt), dtype = complex)
    
    for i in np.arange(int(number_of_iterations)):
        Obs = observable(phi_data[i, :], header, 0)
        phi_numerator[i, :, :] = Obs*phi_tilde[i] #this is the unaveraged numerator
        Obs = observable(chi_data[i, :], header, 1)
        chi_numerator[i, :, :] = Obs*phi_tilde[i]

    #This calculates the expression fully for one initiatlisation
    phi_expectation_observable[n, :, :] = np.mean(phi_numerator, axis = 0)/denominator
    chi_expectation_observable[n, :, :] = np.mean(chi_numerator, axis = 0)/denominator
    
    #error analysis for this file
    for i in np.arange(prot.Nt):
        for j in np.arange(prot.Nt):
            phi_file_error[n, i, j] = jackknife(np.real(phi_numerator[:, i, j]/denominator), jackknife_block_length) + prot.j*jackknife(np.imag(phi_numerator[:, i, j]/denominator), jackknife_block_length)
            chi_file_error[n, i, j] = jackknife(np.real(chi_numerator[:, i, j]/denominator), jackknife_block_length) + prot.j*jackknife(np.imag(chi_numerator[:, i, j]/denominator), jackknife_block_length)
    print("file", n, "completed")

phi_final_expectation = np.mean(phi_expectation_observable, axis = 0)
chi_final_expectation = np.mean(chi_expectation_observable, axis = 0)

phi_final_error = np.sqrt(np.sum(np.real(phi_file_error)**2, axis = 0) + prot.j*np.sum(np.imag(phi_file_error)**2, axis = 0))/np.sqrt(n_files)
chi_final_error = np.sqrt(np.sum(np.real(chi_file_error)**2, axis = 0) + prot.j*np.sum(np.imag(chi_file_error)**2, axis = 0))/np.sqrt(n_files)
location = ""

for i in np.arange(prot.Nt):
    plt.figure(i + prot.Nt)
    plt.errorbar(x_range, np.real(chi_final_expectation[:, i]).flatten(), yerr = np.real(chi_final_error[:, i])/5, label = "real", fmt = '', linestyle = 'None', marker = '.', capsize = 2)
    plt.errorbar(x_range, np.imag(chi_final_expectation[:, i]).flatten(), yerr = np.imag(chi_final_error[:, i])/5, label = "imag", fmt = '', linestyle = 'None', marker = '.', capsize = 2)
    plt.axhline(0,color='black',linewidth = 0.5)
    plt.title('Classical Correlator i = ' + str(x_range[i]))
    plt.xlabel('j')
    plt.ylabel(r'$\langle \chi^{cl}_j \chi^{cl}_i \rangle$')
    plt.legend()
    plt.plot()
    plt.savefig(location + 'classical_chi_'+str(i))
    plt.close()

fig, ax = plt.subplots(prot.Nt, 2, figsize = (5, 10))
fig.tight_layout()
for i in np.arange(prot.Nt):
    ax[i, 0].errorbar(x_range, np.real(phi_final_expectation[:, i]).flatten(), yerr = np.real(phi_final_error[:, i])/5, label = "real", fmt = '', linestyle = 'None', marker = '.', capsize = 5)
    ax[i, 0].errorbar(x_range, np.imag(phi_final_expectation[:, i]).flatten(), yerr = np.imag(phi_final_error[:, i])/5, label = "imag", fmt = '', linestyle = 'None', marker = '.', capsize = 5)
    #ax[i, 0].axhline(0,color='black',linewidth = 0.5)
    #ax[i, 0].set_title('Classical Correlator i = ' + str(x_range[i]))
    #ax[i, 0].set_xlabel('j')
    #ax[i, 0].set_ylabel(r'$\langle \phi^{cl}_j \phi^{cl}_i \rangle$')
    #ax[i, 0].legend()

    ax[i, 1].errorbar(x_range, np.real(chi_final_expectation[:, i]).flatten(), yerr = np.real(chi_final_error[:, i])/5, label = "real", fmt = '', linestyle = 'None', marker = '.', capsize = 5)
    ax[i, 1].errorbar(x_range, np.imag(chi_final_expectation[:, i]).flatten(), yerr = np.imag(chi_final_error[:, i])/5, label = "imag", fmt = '', linestyle = 'None', marker = '.', capsize = 5)
    #ax[i, 1].axhline(0,color='black',linewidth = 0.5)
    #ax[i, 1].set_title('Classical Correlator i = ' + str(x_range[i]))
    #ax[i, 1].set_xlabel('j')
    #ax[i, 1].set_ylabel(r'$\langle \phi^{cl}_j \phi^{cl}_i \rangle$')
    #ax[i, 1].legend()
plt.savefig(location + "subplot")
plt.savefig("single_field_update_test_subplot")