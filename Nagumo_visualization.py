import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from pylab import * 


def nagumo_visualize_2D(grid,v_Jn):
    plt.pcolormesh(grid.spatial_grid.X_J,grid.temporal_grid.t,v_Jn)
    plt.xlabel(r'$x$')
    plt.ylabel(r'$t$')
    plt.colorbar()  
    plt.show()
    
def nagumo_visualize_3D(grid,v_Jn):
    X_J, t = np.meshgrid(grid.spatial_grid.X_J, grid.temporal_grid.t)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(X_J, t, v_Jn, rstride=1, 
                           cstride=1, cmap=cm.jet, linewidth=0.02)          
    
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$t$')
    ax.set_zlabel(r'$v(t,x)$')
    ax.set_zlim(0, 1.5)   
    #fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.show()
    
def nagumo_visualize_potential_change_rate(grid,c):
    plt.plot(grid.temporal_grid.t, c)
    plt.xlabel(r'$t$')
    plt.ylabel(r'temporal rate of change in potential at x = 0')
    plt.show()  
    
def nagumo_visualize_signal_fixed_location(grid,s):
    plt.plot(grid.temporal_grid.t, s)
    plt.xlabel(r'$t$')
    plt.ylabel('potential at x = 0')
    plt.show()   
    
def nagumo_visualize_signal_fixed_time(grid,s):
    plt.plot(grid.spatial_grid.X_J, s) 
    plt.xlabel(r'$x$')
    plt.ylabel(r'solution $v(t_n,x)$')
    plt.show()    

def nagumo_visualize_wave_speed(temp_grid,c_rel,
                                c_sample,c_mean,c_std,M):
    # set up plot
    fig = plt.figure()
     
    # single sample path
    tlim = temp_grid.temp_domain
    
    max1 = np.abs(c_sample - c_rel).max()
    max2 = np.abs(c_mean - c_rel).max()
    max12 = max(max1, max2)
    
    # plot single sample path, c_sample
    ax1 = plt.axes(xlim=tlim, 
                   ylim=(c_rel - 1. * max12, c_rel + 2. * max12))
    line1, = ax1.plot([], [], c='r', 
                      label=r'$sample\, path$')
    line1.set_data(temp_grid.t, c_sample)
    
    # plot c_mean
    ax2 = plt.axes(xlim=tlim, 
                   ylim=(c_rel - 1. * max12, c_rel + 2. * max12))
    
    line2, = ax2.plot([], [], c='k',
                      label=r'$mean\, of\,'+str(M)+'\, samples$')
    line2.set_data(temp_grid.t, c_mean)
    
    plt.fill_between(temp_grid.t, c_mean-2*c_std, c_mean+2*c_std, 
                     color='b', alpha=0.1)
    
    ax2.set_xlabel('$t$')
    ax2.set_ylabel('$C(t)/t$')
    ax2.axhline(c_rel, c='k', ls=':', label=r'$c_{rel}=c_{det}-c_{stoch}$')
    ax2.legend(prop=dict(size=12))
    
    plt.show()
    
def nagumo_visualize_mean_crel(alpha,mean,c_rel,M,xlim,ylim):
    # set up plot
    fig = plt.figure() 
    
    # alpha
    alpha_lim = xlim
    
    # c_rel
    c_rel_lim = ylim
    
    # plot true c_rel   
    ax1 = plt.axes(xlim=alpha_lim, ylim=c_rel_lim)
    line1, = ax1.plot([], [], c='r', 
                      label=r'$c_{rel}(\alpha)$')
    line1.set_data(alpha, c_rel)
    
    # plot mean   
    ax2 = plt.axes(xlim=alpha_lim, ylim=c_rel_lim)
    line2, = ax2.plot([], [], c='k', 
                      label=r'$mean\, of\,'+str(M)+'\, samples$')
    line2.set_data(alpha, mean)

    ax2.set_xlabel('$\alpha$')
    ax2.set_ylabel('$c_{rel}(\alpha), \, mean$')    
    ax2.legend(prop=dict(size=12))
    
    plt.show()
    
    
