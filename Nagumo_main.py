import numpy as np

import Nagumo_model as nm
import Nagumo_helpers as nh
import fcts_semi_implicit_em as nfd

import Nagumo_animation as na
import Nagumo_visualization as nv


if __name__ == "__main__":   
    
    # PARAMETERS
    # model PARAMETERS
    a = 0.1
    b = 1.
    nu = 10.
    alpha = 0.2
       
    # SIMULATION PARAMETERS
    #L = 100
    L = 75.## wave speed plot, 80
    T = 50.
    J = 2000. #J = 50000.## wave speed plot
    #J = 5000
    N = 600. #1000

    # two sided spatial domain
    spatial_two_sided = True

    # boundary condition
    BC = 'n'
    
    ###############################################################
    ################## SETTING INITIAL VALUES #####################
    ###############################################################
    
    # perturbed TW 
    # pert = 0.05; for animation
    pert = 0.05; location_TW = 73.
    pert_TW = nm.NagumoInitialValueTW(pert, location_TW)
    
    # Gaussian pulse
    pert_pulse = 1e-2; location_pulse = 0.5
    pulse = nm.NagumoInitialValuePulse(pert_pulse, location_pulse)
    
    # Kink
    c_1 = 1; c_2 = 5. 
    kink = nm.NagumoInitialValueKink(c_1, c_2)
        
    ##############################################################
    ################### SETTING NAGUMO EQUATION ##################
    ##############################################################
     
    # INITIAL VALUE: perturbed TW 
    spde_nagumo_TW = nm.SPDE_FD_Nagumo(nu,a,b,alpha,pert_TW,BC,
                                       L,J,T,N,spatial_two_sided)
    
    # INITIAL VALUE: Gaussian pulse
    spde_nagumo_pulse = nm.SPDE_FD_Nagumo(nu,a,b,alpha,pulse,BC,
                                         L,J,T,N,spatial_two_sided)

    # INITIAL VALUE: Kink
    spde_nagumo_kink = nm.SPDE_FD_Nagumo(nu,a,b,alpha,kink,BC,
                                         L,J,T,N,spatial_two_sided)
    
    
    ##############################################################
    #################### SIMULATIONS ############################# 
    ##############################################################
    
    # simulation with perturbed TW
    fd_grid, v_Jn = nfd.fd_si_em_simul(spde_nagumo_TW)
        
    #simulation with Gaussian pulse
    #fd_grid, v_Jn = nfd.fd_si_em_simul(spde_nagumo_pulse)
    
    # simulation with kink
    #fd_grid, v_Jn = nfd.fd_si_em_simul(spde_nagumo_kink)
    
    
    ##############################################################
    ################### VISUALIZATION ############################
    ##############################################################    
    
    # visualization in 2D-plot
    #nv.nagumo_visualize_2D(fd_grid, np.transpose(v_Jn))
    
    # visualization in 3D-plot
    #nv.nagumo_visualize_3D(fd_grid, np.transpose(v_Jn))
    
    # visualize wave speed
    #c = nh.nagumo_potential_change_rate(fd_grid, v_Jn)
    #nv.nagumo_visualize_potential_change_rate(fd_grid,c)
    
    # visualize potential
    #s = nh.nagumo_signal_fixed_location(fd_grid, v_Jn)
    #nv.nagumo_visualize_signal_fixed_location(fd_grid, s)
    
    # plot solution
    #s = nh.nagumo_signal_fixed_time(fd_grid, v_Jn, 70)
    #nv.nagumo_visualize_signal_fixed_time(fd_grid, s)
    
    # computing wave speed
    #c_n = nh.nagumo_compute_wave_speed(fd_grid, v_Jn)
    
    # animate wave
    #nagumo_animation = na.NagumoAnimation(spde_nagumo_TW, v_Jn, c_n,
    #                                      [[0,spde_nagumo_TW.c_stoch],[0,1.5]]
    #                                      )
    #nagumo_animation.run_save_animation()
    
