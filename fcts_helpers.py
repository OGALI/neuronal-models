import numpy as np
from scipy import sparse 
import scipy.linalg as LA

from classes_for_grid_handling import TemporalGrid

def composeA(J,BC):
    '''
    ###################### INPUT: #########################
    J       : #of spatial grid points
    BC      : boundary conditions
             'd' - homogeneous Dirichlet; 
             'n' - homogeneous Neumann
    
    ##################### OUTPUT: #########################
    A       : sparse, J-1 x J-1-matrix, 
              if BC=='d', or J x J-matrix, if BC == 'n'
    
    #######################################################
    ''' 
    A = sparse.diags([2.*np.ones((J, )),-np.ones((J-1, )),
                      -np.ones((J-1, ))], [0,-1,1], format='csc')        
    if BC == 'n':  
        A[0,1] = -2.
        A[-1,-2] = -2.
    elif BC == 'd':
        A = A[1:-1,1:-1]
    else:
        raise ValueError("Neither Neumann, nor Dirichlet boundary conditions ")    
    return A    


def generate_bm_sample_path(T,N):
        temp_grid = TemporalGrid(T,N)
        dW = np.random.randn(temp_grid.t.shape[0],)
        dW[0] = 0
        dW = np.sqrt(temp_grid.delta_t) * dW
        return temp_grid, np.cumsum(dW) 


# compute covariance exp(-abs(x)/l)
def compute_exp_cov(h,J,l):
    r = np.arange(J,dtype=float)
    toep = LA.toeplitz(r)
    return (np.exp(-np.abs(toep)*h/l))
    
    
