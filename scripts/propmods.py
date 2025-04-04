import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_heatmap(filename, cur_param, cur_plane, cur_coor):
 # Load data from file
 data = np.loadtxt(filename)
 
 # Extract X, Y coordinates and heatmap values
 x_values = data[:, 0]   # First column → X coordinates
 y_values = x_values   # Second column → Y coordinates
 heatmap = data[:, 1:]   # Remaining columns → Heatmap data
 cur_title = cur_param + " , plane = " + cur_plane + ' , for coor = '+ cur_coor
 
 # Define extent based on X and Y values
 xmin, xmax = x_values.min(), x_values.max()
 ymin, ymax = y_values.min(), y_values.max()
 
 # Plot the heatmap
 plt.imshow(heatmap, interpolation='none', cmap='viridis', extent=[xmin, xmax, ymin, ymax], origin='lower')
 
 # Add labels and colorbar
 plt.colorbar(label="Intensity")
 plt.xlabel("$x[R_*]$")
 plt.ylabel("$y[R_*]")
 plt.title(cur_title)
 
 plt.show()
 



folder = sys.argv[1]

filename_t_xy = "./" + folder + "/propmod_t_xy.dat"
filename_t_yz = "./" + folder + "/propmod_t_yz.dat"
filename_rho_xy = "./" + folder + "/propmod_rho_xy.dat"
filename_rho_yz = "./" + folder + "/propmod_rho_yz.dat"

plot_heatmap(filename_t_xy, 'T', 'xy', '0')
plot_heatmap(filename_t_yz, 'T', 'yz', '0')
plot_heatmap(filename_rho_xy, 'ρ', 'xy', '0')
plot_heatmap(filename_rho_yz, 'ρ', 'yz', '0')

