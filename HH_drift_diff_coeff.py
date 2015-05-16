import numpy as np   

import fcts_helpers as fhelp


#################################################
########## Drift-Coefficients ###################
#################################################   

########### Equation for U #########################
class HH_U_DriftCoeff:
    def __init__(self,param_U,I):
        # PARAMETERS for U
        self.gNa = param_U['gNa']
        self.gK = param_U['gK']
        self.gL = param_U['gL']
        self.ENa = param_U['ENa']
        self.EK = param_U['EK']
        self.EL = param_U['EL']
        
        # EXITATORY SIGNAL
        self.I = I
        
    def eval(self,U,(n,m,h,t,x)):
        f = -(self.gNa*(m**3)*h*(U-self.ENa)
               +self.gK*(n**4)*(U-self.EK)
               +self.gL*(U-self.EL)) + self.I.eval(t,x)  
        return f
    
########### Equations for x = n,m,h ################      
class HH_X_DriftCoeff:
    def __init__(self,param_x,param_gen_env):
        
        self.param_x = param_x
        self.param_gen_env = param_gen_env
        self.TC = param_gen_env['TC']
        self.T_base = param_gen_env['T_base']
        self.Q_10 = param_gen_env['Q_10']
        self.phi = self.Q_10**((self.TC-self.T_base)/10.)
        
        if self.param_x['type'] == 'h':
            self.alpha_x = AlphaH(self.param_x)
            self.beta_x = BetaH(self.param_x)
        else:    
            self.alpha_x = AlphaX(self.param_x)
            self.beta_x = BetaX(self.param_x)    
    
    def eval(self,x,U):
        dx = self.phi*(self.alpha_x.eval(U) * (1. - x)
              -self.beta_x.eval(U) * x)        
        return dx      


#################################################
########## Diffusion-Coefficients ###############
#################################################
class HH_U_DiffusionCoeff:
    '''
    ##################### PURPOSE: ###########################
    Wrapper-class for particular version of additive or 
    multiplicative noise
    
    ##########################################################
    '''
    def __init__(self, noise):
        # instance of AddNoise
        self.noise = noise 


##########################################################
##### HH-diffusion coefficient paper implementation ######
##########################################################    
class HH_X_DiffusionCoeff:
    def __init__(self,sigma_x,J,L,spatial_two_sided):
        self.sigma_x = sigma_x 
        self.J = J
        if spatial_two_sided:
            L = 2. * L
        self.L = L    
        
    def eval(self,x,add_param=None):
        s1 = (self.sigma_x * np.sqrt(self.L) 
               / 24. / np.sqrt(self.J))
        s2 = (self.sigma_x * np.sqrt(self.L) 
               / np.sqrt(2) / 12. / np.sqrt(self.J))     
        
        v1 = x[1:-1] + x[:-2]
        v2 = x[1:-1] - x[:-2]
        v3 = 12. * x[1:-1]**2 - 8. * x[1:-1]**3        
        x1 = (v3 - 3. * v1**2 + v1**3) / v2
        
        v1 = x[1:-1] + x[2:]
        v2 = x[1:-1] - x[2:]
        x2 = (v3 - 3. * v1**2 + v1**3) / v2
        
        x[0] = s2 * ((8.*x[0]**3 - 12.*x[0]**2 + 
                     3.*(x[0] + x[1])**2 - (x[0] + x[1])**3) 
                     / (x[1] - x[0]))  
        x[-1] = s2 * ((-8.*x[-1]**3 + 12.* x[-1]**2 
                - 3.*(x[-2] + x[-1])**2 + (x[-1] + x[-2])**3) 
                      / (x[-1] - x[-2]))
        x[1:-1] = (x1 + x2) * s1
        return x
    
    def eval_with_wp(self,x,wp,add_param=None):
        y = self.eval(x,add_param)
        return y * wp

 
############################################################### 
##### HH-diffusion coefficient standard fd-implementation #####
############################################################### 
class HH_X_DiffusionCoeff_2:
    def __init__(self,sigma_x):
        self.sigma_x = sigma_x
        
    def eval(self,v,add_param=None):
        v[v<0.] = 0.
        v[v>1.] = 0.
        return self.sigma_x*v*(1.-v) 
    
    def eval_with_wp(self,x,wp,add_param=None):
        y = self.eval(x,add_param)
        return y * wp


############################################################### 
### HH-diffusion coefficient with exponential covariance-op.###
############################################################### 
class HH_X_DiffusionCoeff_3:
    def __init__(self,sigma_x,h,J,l):
        self.sigma_x = sigma_x
        self.h  = h
        self.J = J
        self.l = l
        self.B = fhelp.compute_exp_cov(self.h, self.J, self.l)
        
    def eval(self,v,add_param=None):
        v[v<0.] = 0.
        v[v>1.] = 0.
        return self.sigma_x*v*(1.-v) 
    
    def eval_with_wp(self,x,wp,add_param=None):
        y = self.eval(x,add_param)
        return y * np.dot(self.B,wp)
 
 
 
################################################################
########## alpha and beta-coefficient functions in HH ##########
################################################################                
class AlphaX(object):
    def __init__(self,param_x):
        self.ax_1 = param_x['ax_1']
        self.ax_2 = param_x['ax_2']
        self.A_x = param_x['A_x']
        
    def eval(self,U,add_param=None):
        alpha = (self.ax_1 * (U + self.A_x) 
                 / (1. - np.exp(-self.ax_2 * (U + self.A_x))))
        return alpha
    
    
class BetaX(object):
    def __init__(self,param_x):
        self.bx_1 = param_x['bx_1']
        self.bx_2 = param_x['bx_2']
        self.B_x = param_x['B_x']
        
    def eval(self,U,add_param=None):
        beta = self.bx_1 * np.exp(-self.bx_2 * (U + self.B_x))  
        return beta              
        
        
class AlphaH(AlphaX):
    def __init__(self,param_x):
        AlphaX.__init__(self,param_x)
        
    def eval(self,U,add_param=None):
        alpha = (self.ax_1 * np.exp(-self.ax_2*(U+self.A_x)))  
        return alpha
        

class BetaH(BetaX):
    def __init__(self,param_x):
        BetaX.__init__(self,param_x)
        
    def eval(self,U,add_param=None):
        beta = (self.bx_1 /(1.+np.exp(-self.bx_2*(U+self.B_x))))            
        return beta
