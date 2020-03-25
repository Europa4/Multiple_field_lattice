#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 13:44:31 2020

@author: ppxsw1
"""

import numpy as np
import matplotlib.pyplot as plt

import prot

location = []
location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_3/')
location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_4/')
location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder_tau_5/')
location.append('/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/thimble_finder/')

title = []
title.append('3')
title.append('4')
title.append('5')
title.append('15')

for i in np.arange(4):
    data = np.load(location[i] + "data.npz")
    real_data = data['arr_0']
    imag_data = data['arr_1']
    
    
    fig, axs = plt.subplots(3, 9)
    plt.subplots_adjust(left = 0.05, right = 0.95)
    fig.suptitle(r'Thimble with $\tau = $' + title[i])
    for k in np.arange(prot.Nt - 2):
        axs[0, k].hist2d(real_data[k, :], imag_data[k, :])
        axs[0, k].set_title("Site " + str(k + 1))
        axs[0, k].get_xaxis().set_visible(False)
        axs[0, k].get_yaxis().set_visible(False)
        
    axs[1, prot.Nt - 2].hist2d(real_data[prot.Nt - 1, :], imag_data[prot.Nt - 1, :])
    axs[1, prot.Nt - 2].set_title("Site " + str (prot.Nt - 1))
    axs[1, prot.Nt - 2].get_xaxis().set_visible(False)
    axs[1, prot.Nt - 2].get_yaxis().set_visible(False)
    
    for k in np.arange(prot.Nt - 3):
        axs[2, prot.Nt - 3 - k].hist2d(real_data[prot.Nt - 1 + k, :], imag_data[prot.Nt - 1 + k, :])
        axs[2, prot.Nt - 3 - k].set_title("Site " + str(prot.Nt + k))
        axs[2, prot.Nt - 3 - k].get_xaxis().set_visible(False)
        axs[2, prot.Nt - 3 - k].get_yaxis().set_visible(False)
    
    
    axs[2, 0].hist2d(-1*real_data[0,:], -1*imag_data[0,:])
    axs[2, 0].set_title("Site " + str(prot.Nrp + 1))
    axs[2, 0].get_xaxis().set_visible(False)
    axs[2, 0].get_yaxis().set_visible(False)
    
    axs[0, prot.Nt - 2].axis('off')
    for k in np.arange(prot.Nt - 1):
        axs[1, k].axis('off')
    axs[2, prot.Nt - 2].axis('off')

    fig.savefig(location[i] + '_contour_.png')
    