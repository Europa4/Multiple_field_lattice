#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:22:49 2019

@author: ppxsw1
"""

import numpy as np
import pandas as pd
import math as mth
import matplotlib.pyplot as plt
import prot

val = '3'
valdash = '3'
path = '../../Data_vary_3/'
file_name = path + 'phi_0'
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

#Jacobian loading
iJ = np.loadtxt("ijac")
x = iJ[:, 0::2]
y = iJ[:, 1::2]
iJ = x + prot.j*y
for i in np.arange(number_of_iterations):
    phi_data[i, :] = np.real(np.matmul(iJ, phi_data[i, :]))

movement = np.zeros(number_of_iterations - 1)

for i in np.arange(number_of_iterations - 1):
    currant_state = phi_data[i, :]
    next_state = phi_data[i + 1, :]
    state_change = (currant_state - next_state)**2
    movement[i] = np.sqrt(np.sum(state_change))
    
step_size = 500

number_of_averages = int(mth.floor(movement.size/step_size)) + 1
residual = movement.size - step_size*int(mth.floor(movement.size/step_size))
averaged_movement = np.zeros(number_of_averages)

for i in np.arange(number_of_averages - 1):
    averaged_movement[i] = np.mean(movement[step_size*i:step_size*(i + 1)])

averaged_movement[number_of_averages - 1] = np.mean(movement[-1*residual:])

plt.figure(0)
plt.plot(averaged_movement)
plt.xlabel(r'MC step')
plt.ylabel(r'Distance travelled')
plt.title(r'Walking speed for $\delta = $' + val)
plt.savefig(path + 'Walking_speed_'+valdash+'.pdf')

plt.figure(1)
plt.hist(averaged_movement)
plt.xlabel('size of average step for $\delta = $' + val)
plt.ylabel('frequency')
plt.title('Average distance per step')
plt.savefig(path + 'Walking_speed_hist_'+valdash+'.pdf')

plt.figure(2)
tot_square_distance = np.zeros(number_of_iterations - 1)
for i in (np.arange(number_of_iterations - 1) + 1):
    tot_square_distance[i - 1] = np.sum((phi_data[0, :] - phi_data[i, :])**2)
plt.plot(tot_square_distance)
plt.xlabel(r'Monte-Carlo step number')
plt.ylabel(r'Square total distance travelled ($\Delta x^2$)')
plt.title(r'$\delta = $'+val)
plt.savefig(path + 'Total_distance_'+valdash+'.pdf')