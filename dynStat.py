import sys
import numpy as np
import matplotlib.pyplot as plt

plotwidth = 55
plotheight = 40

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


folder = sys.argv[1]
propCfile = "./" + folder + "/dyn_cells.dat"
modCfile = "./" + folder + "/model_grid.dat"
rhoFile = "./" + folder + "/rho.dat"
savefilename = "./" + folder + "/volumes.eps"

outputfolder = sys.argv[1]

dataP = np.genfromtxt(propCfile)
dataRho = np.genfromtxt(rhoFile)

lenProp = file_len(propCfile)
lenModel = file_len(modCfile)
# odpocet prvnich dvou radku
lenModel -= 2

dataMH = np.genfromtxt(modCfile, skip_footer=lenModel)
dataM = np.genfromtxt(modCfile, skip_header=2)

Rstar = dataMH[1]
# seznam obsazenych modelovych bunek
modelcells = []

print(lenModel+1)
propVol = np.zeros(lenModel+1)
modVol = np.zeros(lenModel+1)
radius = dataM[:,1]
print(radius)
# prepocet jednotlivych polomeru na hranice modelovych bunek
boundaries = np.zeros(lenModel)
boundaries[0] = 1.0
for I in range(lenModel-1):
 boundaries[I + 1] = (radius[I+1] + radius[I]) / 2.0
boundaries[lenModel - 1] = radius[lenModel - 1]

print(boundaries)

x=list(range(1,lenProp+2))

# projdeme vsechny propagacni bunky a secteme, jaky objem nalezi jednotlivym modelovym
# bunkam
for I in range(lenProp):
 width = dataP[I,4:7] / Rstar
 modelCell = int(dataP[I, 8])
 if modelCell != lenModel+1:
  modelcells.append(modelCell)
 upCell = dataP[I, 7]
 if upCell == 0:
  curVol = width[0] * width[1] * width[2]
  # print(curVol, modelCell)
  propVol[modelCell - 1] += curVol

# vypocet objemu jednotlivych bunek
for I in range(lenModel-1):
 rPrev = boundaries[I]
 rAct = boundaries[I + 1]
 volume = 4.0/3.0 * np.pi * (rAct**3 - rPrev**3)
 modVol[I] = volume
 print("r=",radius[I], "propvol = ", propVol[I], "modvol = ", volume, " ratio = ",
 propVol[I]/modVol[I])
 if volume == 0:
  print("zero volume for",I)


ratio=propVol/modVol
# print(ratio)
# print(x,propVol)
xmin=1
xmax=10
fntsize=15
# 
fig, ax = plt.subplots(2,figsize=(plotwidth, plotheight))
# 
ax[0].plot(radius[:lenModel], ratio[:lenModel])
ax[0].set_ylabel("$V_{prop}/V_{model}$", fontsize = fntsize)
ax[0].set_xlim(xmin, xmax)
ax[0].tick_params(axis='both', which='major', labelsize=fntsize)
ax[0].yaxis.offsetText.set_fontsize(fntsize)
#ax.hist(modelcells[:lenProp], bins=100)
ax[1].plot(dataRho[:,1], dataRho[:,2],'*')
ax[1].set_xlabel("$R/R_*$", fontsize = fntsize)
ax[1].set_xlim(xmin, xmax)
ax[1].set_ylabel("density", fontsize = fntsize)
ax[1].tick_params(axis='both', which='major', labelsize=fntsize)
ax[1].yaxis.offsetText.set_fontsize(fntsize)
# 
plt.show()
fig.savefig(savefilename)
