import numpy as np


def nagumo_compute_wave_speed(grid,v_Jn): 
    # indices of first appearance of an element >0.5 per column
    space_ind = np.argmax(v_Jn > 0.5, axis = 0)
    x = grid.spatial_grid.X_J[space_ind]
    dx = np.abs(x[1:] - x[:-1])    
    c_n = dx / grid.temporal_grid.delta_t
    return c_n    
    
    
def nagumo_potential_change_rate(grid,v_Jn): 
    # compute column-index corresponding to x=0
    if grid.spatial_grid.spatial_two_sided == True:
        zero_point = np.ceil(grid.spatial_grid.X_J.shape[0] / 2.)
    else:
        zero_point = 0
    
    # compute wave speed
    delta_s = v_Jn[zero_point,1:] - v_Jn[zero_point,:-1]
    c = np.zeros(grid.temporal_grid.t.shape[0]); c[0] = 0.
    c[1:] = grid.temporal_grid.delta_t * delta_s  
    return c         


def nagumo_signal_fixed_location(grid,v_Jn):   
    # compute column-index corresponding to x=0
    if grid.spatial_grid.spatial_two_sided == True:
        zero_point = np.ceil(grid.spatial_grid.X_J.shape[0] / 2.)
    else:
        zero_point = 0              
        
    s = v_Jn[zero_point,:]
    return s

def nagumo_signal_fixed_time(grid,v_Jn,n):
    s = v_Jn[:,n]
    return s
