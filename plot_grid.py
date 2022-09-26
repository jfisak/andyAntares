# a script for plotting a propagation and model grid
import numpy as np
import matplotlib.pyplot as plt
import sys

xmin = -.5e11
xmax = 4.5e11
ymin = xmin
ymax = xmax

plotwidth = 20
plotheight = 20

R_star = 49194829800.e0
R_sun = 6.957E+10

epsilon = 1

in_file = sys.argv[1]
ou_file = sys.argv[2]

input_data = np.genfromtxt(in_file)
nlines = input_data.shape[0]

epsilon2 = input_data[0,3]

model_file = 'model.dat'
model_data = np.genfromtxt(model_file, skip_header=5)
n_mod_lines = model_data.shape[0]


print('nlines = ', nlines)


fig, ax = plt.subplots(1,figsize=(plotwidth, plotheight))
ax.set_aspect('equal', adjustable='box')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

circle = plt.Circle((0,0), R_star, color='blue')
ax.add_patch(circle)

for I in range(nlines):
 x = input_data[I,0] 
 y = input_data[I,1] 
 z = input_data[I,2] 
 wx = input_data[I,3]
 wy = input_data[I,4]
 wz = input_data[I,5]

 if(abs(z) < epsilon):
  x_values = [x, x+wx, x+wx, x, x]
  y_values = [y, y, y+wy, y+wy, y]
  ax.plot(x_values, y_values, color='grey')

for I in range(n_mod_lines):
 x = model_data[I, 0] * R_sun
 y = model_data[I, 1] * R_sun
 z = model_data[I, 2] * R_sun
 if(z < 2*epsilon2 and z > 0):
  ax.scatter(x, y, color='blue', marker=".")
 


plt.show()
fig.savefig(ou_file)





