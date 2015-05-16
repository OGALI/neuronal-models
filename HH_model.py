import numpy as np

import classes_for_spde_handling as spde 
import classes_for_grid_handling as grh

import HH_drift_diff_coeff as hhdd
import fcts_helpers as fhelp
import HH_noise as hhn


class HH_U(spde.SPDE_FD):
    '''
    ###################### PURPOSE: ##########################
    Implements numerical version of the equation for U in 
    the stochastic Hodgkin-Huxley model 
    
    ###################### Attributes: #######################
    '''    
    def __init__(self,
                 param_U,
                 param_disc,
                 BC,
                 add_noise_type,
                 sigma,
                 U_0,
                 I,
                 add_noise_param=None):
                
        # SPDEFiniteDifference PARAMETERS
        f = hhdd.HH_U_DriftCoeff(param_U,I)        
        spde.SPDE_FD.__init__(self,
                              param_U['diff_param'], 
                              f,
                              None,
                              U_0,
                              BC, 
                              param_disc['L'], 
                              param_disc['J'], 
                              param_disc['T'], 
                              param_disc['N'],
                              param_disc['spatial_two_sided'])
        
        # set NOISE
        self.multiplicative_noise = False           
        if add_noise_type == 1:
            self.sigma = hhn.AddNoise1(sigma,
                                self.st_grid.spatial_grid.X_J.shape[0])
        elif add_noise_type == 2:
            Cov = fhelp.compute_exp_cov(self.st_grid.spatial_grid.h,
                                self.st_grid.spatial_grid.X_J.shape[0],
                                add_noise_param['l'])
            self.sigma = hhn.AddNoise2(sigma, Cov)
            print self.sigma.B
        elif add_noise_type == 3:
            pass
        
        # set Wiener process
        self.wp = np.random.randn(self.st_grid.spatial_grid.X_J.shape[0],
                                   self.st_grid.temporal_grid.t.shape[0])
        
class HH_X(spde.SPDE_FD):
    '''
    ###################### PURPOSE: ##########################
    Implements a numerical version of the equations for 
    x=n,m,h in the stochastic Hodgkin-Huxley model 
    
    ##########################################################
    '''    
    def __init__(self,
                 param_x,
                 param_disc,
                 param_gen_env,
                 mult_noise_type,
                 sigma,
                 x_0):
        
        # SPDEFiniteDifference PARAMETERS
        f = hhdd.HH_X_DriftCoeff(param_x,param_gen_env)     
        spde.SPDE_FD.__init__(self,   
                              None,
                              f, 
                              None, # sigma
                              x_0, 
                              None,
                              param_disc['L'],
                              param_disc['J'],
                              param_disc['T'],
                              param_disc['N'],
                              param_disc['spatial_two_sided'])
        
        if mult_noise_type == 1:
            self.sigma = hhdd.HH_X_DiffusionCoeff(sigma,
                                    self.st_grid.spatial_grid.X_J.shape[0],
                                    self.st_grid.spatial_grid.spatial_domain[1],
                                    param_disc['spatial_two_sided'])
        elif mult_noise_type == 2:
            self.sigma = hhdd.HH_X_DiffusionCoeff_2(sigma)
        elif mult_noise_type == 3:
            self.sigma = hhdd.HH_X_DiffusionCoeff_3(sigma,
                                    self.st_grid.spatial_grid.h,
                                    self.st_grid.spatial_grid.X_J.shape[0],
                                    0.05)
        
        # set realization of discrete Wiener process
        self.wp = np.random.randn(self.st_grid.spatial_grid.X_J.shape[0],
                                   self.st_grid.temporal_grid.t.shape[0])
     
     
 
class HH_Coupled:
    '''
    ###################### PURPOSE: ##########################
    Implements numerical version of the coupled system of 
    equations for U and X=(n,m,h) in the stochastic Hodgkin-
    Huxley model 
    
    ##########################################################
    '''    
    # Equation for U
    U = None 
    
    # Equations for n,m,h
    n = None
    m = None
    h = None
    
    # grid
    st_grid = None
      
    def __init__(self,
                 hh_param,
                 add_noise_type,
                 mult_noise_type,
                 U_0,
                 n_0,
                 m_0,
                 h_0,
                 I,
                 add_noise_param=None):    
        
        # set global grid
        pd = hh_param.param_disc
        self.st_grid = grh.SpatioTemporalGrid(pd['L'],
                                              pd['spatial_two_sided'],
                                              pd['J'],
                                              pd['T'],
                                              pd['N'])
        
        # set SPDEs
        # Equation for U
        self.U = HH_U(hh_param.param_U,
                      hh_param.param_disc,  
                      hh_param.BC,
                      add_noise_type,
                      hh_param.param_sigma['sigma_U'],
                      U_0,
                      I,
                      add_noise_param)
                
        # Equations for X=(n,m,h)
        self.n = HH_X(hh_param.param_n,
                      hh_param.param_disc,
                      hh_param.param_gen_env,
                      mult_noise_type,
                      hh_param.param_sigma['sigma_n'],
                      n_0)           

        self.m = HH_X(hh_param.param_m,
                      hh_param.param_disc,
                      hh_param.param_gen_env,
                      mult_noise_type,
                      hh_param.param_sigma['sigma_m'],
                      m_0)
        
        self.h = HH_X(hh_param.param_h,
                      hh_param.param_disc,
                      hh_param.param_gen_env,
                      mult_noise_type,
                      hh_param.param_sigma['sigma_h'],
                      h_0)

                 
            
        
        
        
