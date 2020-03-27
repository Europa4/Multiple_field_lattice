import numpy as np
import pandas as pd
import prot

location = '/run/media/ppxsw1/78fe3857-1897-4617-a65e-83c9aa61be27/long/'
file_name = location + 'phi_' + '1'
data = pd.read_csv(file_name, header = None, skiprows = 1)
header = np.array(np.array(pd.read_csv(file_name, sep='\s+').keys())[0].split(','), dtype = float)
data = np.array(data)
aux_data = data[:, 2*prot.Nrt:]
phi_data = data[:, :2*prot.Nrt]
#This gets it from the CSV state of "real, imag, real, imag" to "real + i*imag, real + i*imag"
temp = np.reshape(phi_data, (phi_data.shape[0], int(phi_data.shape[1]/2), 2))
phi_data = temp[:, :, 0] + prot.j*temp[:, :, 1]

phi_1 = phi_data[0, :(prot.Nt - 1)]
phi_p = np.fft.fft(phi_1)

print("phi_p^d =")
print(np.conj(phi_p))

print("phi_-p =")
print(np.fft.ifft(phi_1))

print("ratio =")
print(np.conj(phi_p)/np.fft.ifft(phi_1))