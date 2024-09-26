# this script is for reading 2D model atmosphere provided by Nicolas Moens
# the input format of the file is npy and radii are provided in format txt
# reading of the file containing radii do not work yet
# input: model file in the npy format
#        text file providing radii
# output: 2D matrix in the txt format
#
# calling: $ $python read_2D_hydronico.py
#
#
import numpy as np

#________________________________________________________
# basic definitions
filename = 'Star98_0051_regridded.npy'
filename_radius = 'Star98_average_prof.txt'
outputfile = 'Star98_hydronico.dat'

# basic indeces
index_r = 0
index_theta = 1

# variable indeces as provided by Nico
index_rho = 0
index_radialv = 7
index_lateralv = 8
index_temp = 11


data = np.load(filename)
data_r = np.loadtxt(filename_radius, delimiter=',')

N_r = data.shape[0]
N_theta = data.shape[1]

# create empty 2D array
output_data = np.zeros((N_r*N_theta, 6))

cur_index = 0
for ind_I in range(N_r): # loop over radii
 for ind_J in range(N_theta): # loop over angles
  output_data[cur_index, 0] = data_r[ind_I, index_r]
  output_data[cur_index, 1] = np.pi/(N_theta - 1) * ind_J
  output_data[cur_index, 2] = data[ind_I, ind_J, index_rho]
  output_data[cur_index, 3] = data[ind_I, ind_J, index_radialv]
  output_data[cur_index, 4] = data[ind_I, ind_J, index_lateralv]
  output_data[cur_index, 5] = data[ind_I, ind_J, index_temp]
  cur_index += 1

np.savetxt(outputfile, output_data, delimiter=' ')
