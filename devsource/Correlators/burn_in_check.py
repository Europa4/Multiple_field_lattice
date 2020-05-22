

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd

import prot

location = "/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/2_free_field_mass_small_step/"
file_number = 0
file_name = location + "phi_" + str(file_number)
data = pd.read_csv(file_name, header = None, skiprows = 1)
header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
data = np.array(data)
phi_data = data[:, :-4]
temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]
chi_data = phi_data[:, prot.Nrt:]
phi_data = phi_data[:, :prot.Nrt]
aux_data = data[:, -4:]
del data
S = aux_data[:, 0]
mct = np.arange(S.size) + 1
plt.figure(0)
plt.plot(mct, S)
plt.xlabel("Monte Carlo time")
plt.ylabel(r"$\Re[S]$")
plt.title("Action stabilisation for iteration " + str(file_number))
plt.savefig(location + "burn_in_" + str(file_number))
plt.savefig("test")
plt.close(0)

plt.figure(1)
phi_data = np.abs(phi_data)
mean_phi_data = np.mean(phi_data, axis = 1)
plt.plot(mean_phi_data)
plt.xlabel("Monte Carlo time")
plt.ylabel("Mean")
plt.title(r"Characteristic value of $\phi$")
plt.savefig("mean_test")
plt.close(1)

std_phi_data = np.std(phi_data, axis = 1)
plt.plot(std_phi_data)
plt.xlabel("Monte Carlo time")
plt.ylabel("STD on the mean")
plt.savefig("std_test")
plt.close(2)