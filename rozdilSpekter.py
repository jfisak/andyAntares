# a python script for the visualisation of differences between two spectra
#
# $ python rozdilSpekter <soubor1> <soubor2> <vystupniSoubor.eps>
#
# soubor je ve tvaru
# vlnova delka  tok
# Angstroms     Erg/cm^2/...
# 


import numpy as np
import pylab
import matplotlib.pyplot as plt
import sys

plotwidth = 55
plotheight = 40

# axes limits
xmin0 = 800
xmax0 = 1500
xmin1= 4200
xmax1 = 5000
ymin1= 0e1
ymax1 = 1.5e39
xmin2 = 6250
xmax2 = 6750
ymin2= 0e1
ymax2 = 0.7e39

savefilename='tc_rozdil.eps'


spec1= input('First datafile name:')
spec2= input('Second datafile name:')
if spec1 == '':
 spec1 = 'spec1'
if spec2 == '':
 spec1 = 'spec2'

def find_index(spektrum, lambdaw, ncolumns, n_rows):
 for J in range(n_rows):
  if lambdaw >= spektrum[J,0]:
   hledanyIndex = J
   break
  hledanyIndex = -1
 return hledanyIndex

def lininterpolation(spektrum, lambdaw, ncolumns, n_rows):
 for J in range(n_rows):
  if lambdaw >= spektrum[J,0]:
   hledanyIndex = J
   break

 if hledanyIndex == 1:
  hledanyIndex = 2
 
 I0 = hledanyIndex -1
 I1 = hledanyIndex

 
 wale0 = spektrum[I0,0]
 wale1 = spektrum[I1,0]
 func0 = spektrum[I0,1]
 func1 = spektrum[I1,1]

 a = (func1 - func0) / (wale1 - wale0)
 b = (func0 * wale1 - func1 * wale0) / (wale1 - wale0)

 funkcniHodnota = a * lambdaw + b

 return funkcniHodnota

#__________________________________________________________________________________
# nacist nazvy prislusnych souboru
# jsou to parametry pri volani skriptu

soubor1 = sys.argv[1]
soubor2 = sys.argv[2]
outputfile = sys.argv[3]


# načtení dat z jednotlivých souborů
soubor1data = np.genfromtxt(soubor1)
soubor2data = np.genfromtxt(soubor2)

n_one = soubor1data.shape[0]
n_two = soubor2data.shape[0]
ncolumns = soubor1data.shape[1]

# nastavení počátečního a konečného indexu
if soubor2data[0,0] > soubor1data[0,0]:
 end_lambda = soubor1data[0,0]
 I_1_init = 0
 I_2_init = find_index(soubor2data, end_lambda, ncolumns, n_two)
else:
 end_lambda = soubor2data[0,0]
 I_1_init = find_index(soubor1data, end_lambda, ncolumns, n_one)
 I_2_init = 0

if soubor2data[n_two-1, 0] > soubor1data[n_one-1, 0]:
 init_lambda = soubor2data[n_two-1, 0]
 I_2_end = n_two-1
 I_1_end = find_index(soubor1data, init_lambda, ncolumns, n_one)
else:
 init_lambda = soubor1data[n_one-1, 0]
 I_2_end = find_index(soubor2data, init_lambda, ncolumns, n_two)
 I_1_end = n_one-1


# print('I_2_init = ', I_2_init, ' I_2_end = ', I_2_end)
# print('I_1_init = ', I_1_init, ' I_1_end = ', I_1_end)

n_rows_output = I_1_end - I_1_init + 1
diffdata = np.zeros((n_rows_output, 4))

for I in range (I_1_init, I_1_end):
 wavelength = soubor1data[I,0]
 flux1 = soubor1data[I,1]
 flux2 = lininterpolation(soubor2data, wavelength, ncolumns, n_two)
 # print('flux1 = ', flux1, ' flux2 = ', flux2)
 if flux1 > 0 and flux2 > 0:
  diffdata[I,0] = wavelength
  diffdata[I,1] = flux1
  diffdata[I,2] = flux2
  diffdata[I,3] = (flux2-flux1)/flux2
 else:
  diffdata[I,0] = wavelength
  diffdata[I,1] = flux1
  diffdata[I,2] = flux2
  diffdata[I,3] = 0
  # print('I = ', I, ' flux2 = ', diffdata[I,2])

#__________________________________________________________________________________
# plot of spectra and differences

fig, ax = plt.subplots(ncols=3, nrows=2, figsize=(plotwidth, plotheight))
ax[0,0].plot(diffdata[:,0], diffdata[:,1], linewidth=2.0, label=spec1)
ax[0,0].plot(diffdata[:,0], diffdata[:,2], linewidth=3.0, label=spec2)
ax[0,0].set_xlim(xmin0, xmax0)
ax[0,0].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[0,0].yaxis.offsetText.set_fontsize(60)
ax[1,0].set_xlabel("$\lambda [\AA ]$", fontsize = 60)
ax[0,0].set_ylabel("flux/erg/s/$\AA$", fontsize = 60)

ax[1,0].set_xlim(xmin0, xmax0)
ax[1,0].plot(diffdata[:,0], diffdata[:,3], linewidth=.5, color="black")
ax[1,0].set_ylabel("$(I_2-I_1)/I_2$", fontsize = 60)
ax[1,0].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[1,0].yaxis.offsetText.set_fontsize(60)
ax[1,0].set_ylim(-1, 1)
# ax2.set_aspect("equal")
# ax2.grid(which='major', axis='y')
ax[1,0].axhline(y=0, color='k')

####

ax[0,1].plot(diffdata[:,0], diffdata[:,1], linewidth=2.0, label=spec1)
ax[0,1].plot(diffdata[:,0], diffdata[:,2], linewidth=3.0, label=spec2)
ax[0,1].set_xlim(xmin1, xmax1)
ax[0,1].set_ylim(ymin1, ymax1)
ax[0,1].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[0,1].yaxis.offsetText.set_fontsize(60)
ax[1,1].set_xlabel("$\lambda [\AA ]$", fontsize = 60)

#ax2=ax.twinx()
ax[1,1].set_xlim(xmin1, xmax1)
ax[1,1].plot(diffdata[:,0], diffdata[:,3], linewidth=.5, color="black")
ax[1,1].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[1,1].yaxis.offsetText.set_fontsize(60)
ax[1,1].set_ylim(-1, 1)
# ax2.set_aspect("equal")
# ax2.grid(which='major', axis='y')
ax[1,1].axhline(y=0, color='k')

####

ax[0,2].plot(diffdata[:,0], diffdata[:,1], linewidth=2.0, label=spec1)
ax[0,2].plot(diffdata[:,0], diffdata[:,2], linewidth=3.0, label=spec2)
ax[0,2].set_xlim(xmin2, xmax2)
ax[0,2].set_ylim(ymin2, ymax2)
ax[0,2].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[0,2].yaxis.offsetText.set_fontsize(60)
ax[0,2].legend(loc=1, prop={'size': 60})
ax[1,2].set_xlabel("$\lambda [\AA ]$", fontsize = 60)

#ax2=ax.twinx()
ax[1,2].set_xlim(xmin2, xmax2)
ax[1,2].plot(diffdata[:,0], diffdata[:,3], linewidth=.5, color="black")
ax[1,2].tick_params(axis='both', which='major', labelsize=60, width=2, length=10)
ax[1,2].yaxis.offsetText.set_fontsize(60)
ax[1,2].set_ylim(-1, 1)
# ax2.set_aspect("equal")
# ax2.grid(which='major', axis='y')
ax[1,2].axhline(y=0, color='k')
fig.show()
fig.savefig(outputfile)
