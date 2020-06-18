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

outputfolder = sys.argv[1]

dataP = np.genfromtxt(propCfile)

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
x=list(range(1,lenProp+2))

# projdeme vsechny propagacni bunky a secteme, jaky objem nalezi jednotlivym modelovym
# bunkam
for I in range(lenProp):
 width = dataP[I,4:7] / Rstar
 modelCell = int(dataP[I, 8])
 if modelCell != 902:
  modelcells.append(modelCell)
 upCell = dataP[I, 7]
 if upCell == 0:
  curVol = width[0] * width[1] * width[2]
  # print(curVol, modelCell)
  propVol[modelCell - 1] += curVol

# vypocet objemu jednotlivych bunek
for I in range(lenModel-1):
 if I == 0:
  rPrev = 1.E0
 else:
  rPrev = dataM[I-1,1]
 rAct = dataM[I,1]
 volume = 4.0/3.0 * np.pi * (rAct**3 - rPrev**3)
 modVol[I] = volume
 if volume == 0:
  print("zero volume for",I)

ratio=propVol/modVol
# print(ratio)
# print(x,propVol)
# 
fig, ax = plt.subplots(figsize=(plotwidth, plotheight))
# 
ax.plot(x[:lenModel], ratio[:lenModel])
# ax.hist(modelcells[:lenProp], bins=100)
# 
plt.show()
