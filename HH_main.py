import numpy as np

import HH_animation as hha
import HH_parameters as hhp 
import HH_initial_values as hhinit    
import HH_visualization as hhv
import HH_model as hhm
import HH_fd2 as hhfd


if __name__=="__main__":

    ############ DISCRETIZATION parameters ##############
    L = 10.
    T = 20. # duration of action potential (ms) 
    J = 500.
    N = 500
        
    # two sided spatial domain
    spatial_two_sided = True
    
    param_disc = {
                  'L':L,
                  'T':T,
                  'J':J,
                  'N':N,
                  'spatial_two_sided':spatial_two_sided 
                  }    

    # boundary condition
    BC = 'n'
    
    ############## NOISE parameters ######################
    sigma_U = 0.
    sigma_n = 0.
    sigma_m = 0.
    sigma_h = 0.
    
    param_sigma = {
                   'sigma_U':sigma_U,
                   'sigma_n':sigma_n,
                   'sigma_m':sigma_m,
                   'sigma_h':sigma_h
                   }
    
    # set HHData object
    hh_param = hhp.HHData(hhp.param_gen_env,
                          hhp.param_U,
                          hhp.param_n,
                          hhp.param_m,
                          hhp.param_h,
                          param_sigma,
                          param_disc,
                          BC)
    
    ################################################################
    ################## SETTING INITIAL VALUES ######################
    ################################################################    
    rs = hhinit.RestingStates(hh_param)
    rp = rs.compute_resting_U()
    (n_inf,m_inf,h_inf) = rs.compute_resting_X(rp)
        
    U_0 = hhinit.HHInitialEquilibrium(rp)
    n_0 = hhinit.HHInitialEquilibrium(n_inf)
    m_0 = hhinit.HHInitialEquilibrium(m_inf)
    h_0 = hhinit.HHInitialEquilibrium(h_inf)
    
    ################################################################
    ################## SETTING EXCITATORY SIGNAL ###################
    ################################################################
    pert1 = 50.
    pert2 = .5
    loc = L
    
    # time duration of excitatory shock (ms) 
    signal_start = 0.2*T
    signal_end = 0.6*T
    
    I = hhinit.ExcitatorySignal(signal_start,
                                signal_end,
                                pert1,
                                pert2,
                                loc)
    
    ################################################################
    ################## SETTING ADDITIVE NOISE ######################
    ################################################################
    
    # ADDITIVE NOISE for U
    add_noise_type = 1
    
    add_noise_param = {
                       'l':.01
                       }
    
    ################################################################
    ############### SETTING MULTIPLICATIVE NOISE ###################
    ################################################################
    mult_noise_type = 2
    
    
    ################################################################
    ################## SETTING HH EQUATION #########################
    ################################################################
    spde_hh = hhm.HH_Coupled(hh_param,
                             add_noise_type, 
                             mult_noise_type,
                             U_0, 
                             n_0, 
                             m_0, 
                             h_0, 
                             I,
                             add_noise_param)
    
    ################################################################
    ################## SIMULATIONS #################################
    ################################################################
    theta = 1.
    
    st_grid,U_Jn,n_Jn,m_Jn,h_Jn = hhfd.fd_simul_hodgkin_huxley(spde_hh,
                                                               theta)
                                                        
    ################################################################
    ################## VISUALIZATION ###############################
    ################################################################    
    hhv.HH_visualize_2D(st_grid, np.transpose(U_Jn))
    #hhv.HH_visualize_signal_fixed_time(st_grid, U_Jn[:,50])
    #hhv.HH_visualize_signal_fixed_location(st_grid, U_Jn[0,:])
    
    hhv.HH_visualize_2D(st_grid, np.transpose(n_Jn))
    hhv.HH_visualize_2D(st_grid, np.transpose(m_Jn))
    hhv.HH_visualize_2D(st_grid, np.transpose(h_Jn))
    
    ################################################################
    ################## ANIMATION ###################################
    ################################################################    
    plot_ranges = [[0.,1.],[-100.,70.]]
    
    hh_animation = hha.HHAnimation(spde_hh, 
                                   U_Jn, 
                                   n_Jn,
                                   m_Jn,
                                   h_Jn,
                                   plot_ranges)
    hh_animation.run_save_animation() 
    
