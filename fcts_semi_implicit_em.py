import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as splinalg

  
def compose_iteration_matrix_fd(A,J,diff_param,delta_t,h):
    # iteration matrix
    EE = sparse.eye(J) + diff_param * delta_t * A/h/h 
    return EE
    

def fd_si_em_single_step(spde,v_n,xi,EE,
                         add_param_f=None,
                         add_param_sigma=None):
    '''
    #################### PURPOSE: ################################
    Performs one step in a semi-implicit Euler-Maruyama 
    finite-difference scheme
       
    #################### INPUT: ##################################
    spde         : instance of a class implementing a numerical 
                   SPDE 
    v_n          : current iteration vector
    xi           : random sample
    EE           : Iteration matrix
    add_param_f  : length-N tuple containing all other relevant 
                   data for evaluating coefficient function f in 
                   spde
    add_param_sigma  : length-M tuple containing all relevant
                       data for evaluating coefficient function
                       sigma in spde              
    
    #################### OUTPUT: #################################
    v_n        : subsequent iteration vector
    
    ##############################################################
    '''    
    delta_t = spde.st_grid.temporal_grid.delta_t
    f_n = spde.f.eval(v_n,add_param_f)
    if spde.multiplicative_noise:
        sigma_n_W = spde.sigma.eval_with_wp(v_n,xi,add_param_sigma)
    else:
        sigma_n_W = spde.sigma.eval_with_wp(xi)
    v_n = (splinalg.spsolve(EE, v_n + delta_t * f_n
                            + np.sqrt(delta_t) * sigma_n_W))
    return v_n 


def fd_si_em_simul(spde):
    '''
    #################### PURPOSE: ################################
    Generates one sample path for a SPDE using finite differences 
    (in space) and a semi-implicit Euler-Maruyama scheme (in time)
    
    #################### INPUT: ##################################
    spde      : instance of a class implementing a numerical SPDE
     
    #################### OUTPUT: #################################
    st_grid   : instance of class SpatioTemporalGrid
    v_Jn      : matrix, columns approximate one sample
                [v(t_n,x_0),...,v(t_n,x_J)]^T
    
    ##############################################################
    '''

    # spatial and temporal grid
    st_grid = spde.st_grid
    temporal_grid = st_grid.temporal_grid
    spatial_grid = st_grid.spatial_grid
     
    l = temporal_grid.t.shape[0]
    J = spatial_grid.X_J.shape[0]

    delta_t = temporal_grid.delta_t
    h = spatial_grid.h
    
    # set parameters 
    A = spde.A
    diff_param = spde.diff_param 
    BC = spde.BC
    
    # Wiener process
    xi = spde.wp
    
    # compose iteration matrix
    EE = compose_iteration_matrix_fd(A,J,diff_param,delta_t,h,BC)
       
    # initialize 
    v_Jn = np.zeros((J,l)) 
    v_Jn[:,0] = spde.init_value.eval(spatial_grid.X_J)                
    for n in xrange(l-1):
        v_Jn[:,n+1] = fd_si_em_single_step(spde, 
                                           v_Jn[:,n],  
                                           xi[n], EE)                       
    return st_grid, v_Jn


    
