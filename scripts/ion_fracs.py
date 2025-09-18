import numpy as np
import matplotlib.pyplot as plt
import sys

plotwidth = 20
plotheight = 20

in_folder = sys.argv[1]

in_file = in_folder + '/ionFracs.dat'
temp_file = in_folder + '/tempStruct.dat'

def find_element(atom_number):
 if(atom_number == 1):
  elname = 'H'
 elif(atom_number == 2):
  elname = 'He'
 elif(atom_number == 3):
  elname = 'Li'
 elif(atom_number == 4):
  elname = 'Be'
 elif(atom_number == 5):
  elname = 'B'
 elif(atom_number == 6):
  elname = 'C'
 elif(atom_number == 7):
  elname = 'N'
 elif(atom_number == 8):
  elname = 'O'
 else:
  elname = 'bablbam'
 return elname

def find_ion(ion_number):
 if(ion_number == 1):
  ionname = 'I'
 elif(ion_number == 2):
  ionname = 'II'
 elif(ion_number == 3):
  ionname = 'III'
 elif(ion_number == 4):
  ionname = 'IV'
 elif(ion_number == 5):
  ionname = 'V'
 elif(ion_number == 6):
  ionname = 'VI'
 elif(ion_number == 7):
  ionname = 'VII'
 elif(ion_number == 8):
  ionname = 'VIII'
 else:
  ionname = '666'
 return ionname


ion_data = np.genfromtxt(in_file)
t_data = np.genfromtxt(temp_file)
nrows = ion_data.shape[0]
elements = []
n_modgrid_cells = ion_data[nrows - 1, 0]

# numbe of elements
n_elements = 0
last_element = 0
first_modindex = ion_data[0,0]
for I in range (nrows):
 mod_index = ion_data[I,0]
 if(mod_index == first_modindex):
  cur_element = ion_data[I,1]
  if(cur_element != last_element):
   elements.append(cur_element)
   last_element = cur_element
   n_elements += 1
 if(mod_index != first_modindex):
  break

cur_ion = 0
print(n_elements)
for I in range(n_elements):
 element = elements[I]
 n_ions = element + 1
 fractions = np.zeros((int(n_ions), int(n_modgrid_cells)))
 radii = np.zeros(int(n_modgrid_cells))
 for J in range(nrows):
  cur_mgi = ion_data[J, 0]
  cur_element = int(ion_data[J, 1])
  if(cur_element == element):
   cur_ion = ion_data[J, 2]
   if(np.isnan(cur_ion)):
    continue
   cur_frac = ion_data[J, 3]
   cur_radius = ion_data[J, 4]
   fractions[int(cur_ion) - 1, int(cur_mgi) - 1] = cur_frac
   radii[int(cur_mgi) - 1] = cur_radius

 fig, ax = plt.subplots(1,figsize=(plotwidth, plotheight))
 ax2 = ax.twinx()
 ax.set_yscale('log')
 element_name = find_element(element)
 print(element_name)

 for I in range(int(n_ions)):
  ion_name = find_ion(I+1)
  full_ion_name = element_name + ' ' + ion_name
  ax.plot(radii[:], fractions[I, :], 'o', label=full_ion_name)
 ax2.plot(t_data[:,0], t_data[:,1], label='$T$')
 ax.legend(loc=1)
 plt.show()



 del fractions

