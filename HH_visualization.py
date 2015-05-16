import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import * 

def HH_visualize_2D(grid,v_Jn):
    plt.pcolormesh(grid.spatial_grid.X_J,grid.temporal_grid.t,v_Jn)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$t$')
    plt.colorbar()  
    plt.show()
    
def HH_visualize_signal_fixed_time(grid,s):
    plt.plot(grid.spatial_grid.X_J, s) 
    plt.xlabel(r'$x$')
    plt.ylabel(r'solution $v(t_n,x)$')
    plt.show()     
    
def HH_visualize_signal_fixed_location(grid,s):
    plt.plot(grid.temporal_grid.t, s)
    plt.xlabel(r'$t$')
    plt.ylabel('potential at x = 0')
    plt.show()    
