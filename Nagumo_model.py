import numpy as np
import classes_for_spde_handling as spde


class SPDE_FD_Nagumo(spde.SPDE_FD): 
    '''
    #################### PURPOSE: ##########################
    Sub-class of class SPDE_FD for finite differences 
    approximation of stochastic Nagumo's equation
    
    #################### ATTRIBUTES: #######################
    '''     
    # REACTION TERM: parameters
    a = None
    b = None
    
    # NOISE AMPLITUDE
    alpha = None
    
    # DIFFUSION: parameters
    nu = None
    nu_bar = None
    
    # WAVE SPEEDS
    c_det = None
    c_stoch = None
     
    def __init__(self,nu,a,b,alpha,v_0,BC,
                 L,J,T,N,spatial_two_sided):
        
        # REACTION TERM: parameters
        self.a = a 
        self.b = b 
        
        # NOISE AMPLITUDE
        self.alpha = alpha
        
        # DIFFUSION: parameters
        self.nu = nu 
        self.nu_bar = (nu*(1+b*alpha**2))        
       
        # WAVE SPEEDS
        # deterministic case
        self.c_det = (np.sqrt(2.*self.b*self.nu_bar) 
                      *(0.5-self.a))
         
        # stochastic case
        self.c_stoch = (np.sqrt(2.*self.b*self.nu) 
                        * (0.5-self.a))       
       
        # Base-class: SPDE_FD
        # PARAMETERS
        f = NagumoDriftCoeff(self.a, self.b)
        sigma = NagumoDiffusionCoeff(self.alpha * self.b)
        
        spde.SPDE_FD.__init__(self, self.nu_bar,
                              f, sigma, v_0, BC,
                              L, J, T, N,
                              spatial_two_sided)

        # Wiener process
        self.wp = np.random.randn(self.st_grid.temporal_grid.t.shape[0],)


######################################################             
########### Drift-/Diffusion-Coefficients ############
######################################################        
class NagumoDriftCoeff:  
    def __init__(self,a,b):
        self.a = a
        self.b = b
    
    def eval(self,v,add_param=None):
        return self.b*v*(1.-v)*(v-self.a)

    
class NagumoDiffusionCoeff:       
    def __init__(self,sigma_0):
        self.sigma_0 = sigma_0
        
    def eval(self,v,add_param=None):
        v[v<0.] = 0.
        v[v>1.] = 0.
        return self.sigma_0*v*(1.-v)    
    
    def eval_with_wp(self,v,wp,add_param=None):
        return self.eval(v,add_param) * wp


######################################################             
################ Initial Values ######################
######################################################          
class NagumoInitialValueTW:  
    def __init__(self,perturbation,location):
        self.perturbation = perturbation
        self.location = location 
                
    def eval(self,v,add_param=None):
        return 1. /(1. + np.exp((-self.perturbation 
                                 * (v - self.location) / np.sqrt(2))))
 
    
class NagumoInitialValuePulse:   
    def __init__(self,perturbation,location):
        self.perturbation = perturbation
        self.location = location
            
    def eval(self,v,add_param=None):
        return np.exp(-(self.perturbation * (v - self.location)**2.))    
    
    
class NagumoInitialValueKink:   
    def __init__(self,c_1,c_2):
        self.c_1 = c_1
        self.c_2 = c_2
        
    def eval(self,v,add_param=None):
        a = -np.ones((v.shape[0],))
        a[np.abs(v) <= self.c_2] = -self.c_1
        a[np.abs(v) > self.c_2] = self.c_1        
        a[np.abs(v) > 3. * self.c_2] = 0. 
        return a
        #return self.c_1 * np.sin(v)
