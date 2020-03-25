#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 12:10:21 2020

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

import prot

n_files = 200
x_range = np.arange(prot.Nt - 2) + 1
location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/4_site/'
real_data = np.zeros((prot.Nrt, 0), dtype = float)
imag_data = np.copy(real_data)
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
    real_data = np.concatenate((real_data, np.real(np.transpose(phi_data))), axis=1)
    imag_data = np.concatenate((imag_data, np.imag(np.transpose(phi_data))), axis=1)

#for i in np.arange(prot.Nt):
#    plt.figure(i)
#    plt.plot(real_data[i, :], imag_data[i, :], ',')
#    plt.show()
    
np.savez(location + "data", real_data, imag_data)

for i in np.arange(prot.Nrt):
    plt.figure(i)
    #plt.hist2d(real_data[i, :], imag_data[i, :], bins=20, range=[[-3, 3], [-3, 3]])
    plt.hist2d(real_data[i,:], imag_data[i,:], bins = 100)
    plt.colorbar()
    plt.xlabel(r'$real(\phi)$')
    plt.ylabel(r'$imag(\phi)$')
    plt.title('thimble at site ' + str(i))
    plt.savefig(location + '_scatter_' + str(i))   




#n_files = 8
#x_range = np.arange(prot.Nt - 2) + 1
#location = []
#location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_3/')
#location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_4/')
#location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_5/')
#location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder/')
#real_data = np.zeros((prot.Nrt, 0), dtype = float)
#imag_data = np.copy(real_data)
#for k in np.arange(4):
#    for n in np.arange(n_files):
#        #data location
#        file_name = location[k] + 'phi_' + str(n)
#        #This is all just using pandas to load the CSV files, it's completely compatable with the outputs from the C++ code
#        #It outputs it as a triplet of numpy arrays, header, aux_data and phi_data
#        data = pd.read_csv(file_name, header = None, skiprows = 1)
#        header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
#        data = np.array(data)
#        aux_data = data[:, 2*prot.Nrp:]
#        phi_data = data[:, :2*prot.Nrp]
#        #This gets it from the CSV state of "real, imag, real, imag" to "real + i*imag, real + i*imag"
#        temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
#        phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]
#        number_of_iterations = phi_data.shape[0]
#        real_data = np.concatenate((real_data, np.real(np.transpose(phi_data))), axis=1)
#        imag_data = np.concatenate((imag_data, np.imag(np.transpose(phi_data))), axis=1)
#        
#        for i in np.arange(prot.Nrt):
#            plt.figure(i)
#            plt.subplot(2, 2, k + 1)
#            plt.hist2d(real_data[i,:], imag_data[i,:])
#            #plt.colorbar()
#            plt.xlabel(r'$real(\phi)$')
#            plt.ylabel(r'$imag(\phi)$')
#            
#            
#        real_data = np.zeros((prot.Nrt, 0), dtype = float)
#        imag_data = np.copy(real_data)
            
        
    
#    for i in np.arange(prot.Nrt):
#        plt.figure(i)
#        #plt.hist2d(real_data[i, :], imag_data[i, :], bins=20, range=[[-3, 3], [-3, 3]])
#        plt.hist2d(real_data[i,:], imag_data[i,:])
#        plt.colorbar()
#        plt.xlabel(r'$real(\phi)$')
#        plt.ylabel(r'$imag(\phi)$')
#        plt.title('thimble at site ' + str(i))
#        plt.savefig(location + '_scatter_' + str(i))   