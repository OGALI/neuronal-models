import numpy as np

import fcts_semi_implicit_em as si_em
import fcts_theta_em as tem


def fd_simul_hodgkin_huxley(spde_hh,theta):    
    # Equation for U
    U = spde_hh.U
    
    # Equations for n,m,h
    n = spde_hh.n
    m = spde_hh.m
    h = spde_hh.h
    
    # spatial and temporal grid
    st_grid = spde_hh.st_grid
    temporal_grid = st_grid.temporal_grid
    spatial_grid = st_grid.spatial_grid
    
    l = temporal_grid.t.shape[0]
    J = spatial_grid.X_J.shape[0]
    
    delta_t = temporal_grid.delta_t
    delta_h = spatial_grid.h
    
    # parameters for U
    A = U.A
    diff_param = U.diff_param
        
    # Wiener processes
    xi_U = U.wp
    xi_n = n.wp
    xi_m = m.wp
    xi_h = h.wp
    
    # compose iteration matrix
    EE = si_em.compose_iteration_matrix_fd(A, J, diff_param, 
                                           delta_t, delta_h)
    
    # initialize
    U_Jn = np.zeros((J,l))
    U_Jn[:,0] = U.init_value.eval(spatial_grid.X_J)
    #print U_Jn[:,0]
    
    n_Jn = np.zeros((J,l))
    m_Jn = np.zeros((J,l))
    h_Jn = np.zeros((J,l)) 
    n_Jn[:,0] = n.init_value.eval(spatial_grid.X_J)
    m_Jn[:,0] = m.init_value.eval(spatial_grid.X_J)
    h_Jn[:,0] = h.init_value.eval(spatial_grid.X_J)
    
    for i in xrange(l-1):
        print i
        # first: one step for U
        U_Jn[:,i+1] = si_em.fd_si_em_single_step(U, U_Jn[:,i], 
                                xi_U[:,i], EE, 
                                add_param_f=(n_Jn[:,i],m_Jn[:,i],h_Jn[:,i],
                                            temporal_grid.t[i],
                                            spatial_grid.X_J))
        # second: one step for each n,m,h in X
        #right hand sides    
        rs_n = (n_Jn[:,i]+(1.-theta)*delta_t*n.f.eval(n_Jn[:,i],U_Jn[:,i])
            +np.sqrt(delta_t)*n.sigma.eval_with_wp(n_Jn[:,i],xi_n[:,i]))
 
        rs_m = (m_Jn[:,i]+(1.-theta)*delta_t*m.f.eval(m_Jn[:,i],U_Jn[:,i])
            +np.sqrt(delta_t)*m.sigma.eval_with_wp(m_Jn[:,i],xi_m[:,i]))

        rs_h = (h_Jn[:,i] + (1.-theta)*delta_t*h.f.eval(h_Jn[:,i],U_Jn[:,i])
            +np.sqrt(delta_t)*h.sigma.eval_with_wp(h_Jn[:,i],xi_h[:,i]))
        
        n_Jn[:,i+1] = tem.theta_em_single_step(n, n_Jn[:,i], 
                          (theta,delta_t,rs_n,U_Jn[:,i+1]))
        m_Jn[:,i+1] = tem.theta_em_single_step(m, m_Jn[:,i], 
                          (theta,delta_t,rs_m,U_Jn[:,i+1]))
        h_Jn[:,i+1] = tem.theta_em_single_step(h, h_Jn[:,i],
                          (theta,delta_t,rs_h,U_Jn[:,i+1]))
                        
    return st_grid, U_Jn, n_Jn, m_Jn, h_Jn 
