import classes_for_grid_handling as grh
from fcts_helpers import composeA

class SPDE(object):
    '''
    ################# PURPOSE: #########################
    Base class for numerical versions of stochastic 
    reaction diffusion equations 
        
    ################# ATTRIBUTES: ######################
    '''
    # DIFFUSION
    # diffusion parameter
    diff_param = None 
    
    # boundary condition
    BC = None
                
    # discrete Laplacian
    A = None
        
    # REACTION TERM 
    f = None
        
    # NOISE
    multiplicative_noise = True
    
    # stochastic diffusion coefficient
    sigma = None
        
    # Wiener process
    wp = None
        
    # INITIAL VALUE
    init_value = None
    
    def __init__(self,diff_param,f,sigma,
                 init_value,A=None,BC=None):
        # set diffusion
        self.diff_param = diff_param
        self.BC = BC        
        self.A = A
        
        # set reaction term
        self.f = f 
        
        # set noise
        self.sigma = sigma
        
        # initial value
        self.init_value = init_value
    
 
class SPDE_FD(SPDE): 
    '''
    #################### PURPOSE: #######################
    Sub-class of class SPDE implementing a finite 
    difference approximation
    
    #################### ATTRIBUTES: ####################
    ''' 
    def __init__(self,diff_param,f,sigma,init_value,BC,
                 L,J,T,N,spatial_two_sided):  
        # GRID
        self.st_grid = grh.SpatioTemporalGrid(L,
                                              spatial_two_sided,
                                              J,T,N)
        # LAPLACIAN 
        # with boundary conditions BC
        if diff_param != None:
            A = composeA(self.st_grid.spatial_grid.X_J.shape[0],BC)
        else:
            A = None
                
        # set general attributes of SPDE
        SPDE.__init__(self,diff_param,f,sigma,init_value,A,BC)


##############################################################
####### Base-Classes for Drift/Diffusion/Noise wrappers ######
############################################################## 
class DriftCoeff:
    def __init__(self,drift):
        self.drift = drift 
        
        
class DiffusionCoeff:
    def __init__(self,diffusion):
        self.diffusion = diffusion        


class Noise(object):
    '''
    ##################### PURPOSE: ###########################
    Common Base-Class for AddNoise and MultNoise classes
    
    ##########################################################    
    '''  
    def __init__(self):
        pass 
    
    def eval_with_wp(self):
        pass 
          
        
class AddNoise(Noise):
    '''
    ##################### PURPOSE: ###########################
    Base-class for additive noise; implements an efficient 
    evaluation of additive noise with a realization of a 
    Wiener process wp
    
    ##########################################################
    '''
    def __init__(self,sigma,B): 
        self.sigma = sigma
        self.B = B
        
    def eval_with_wp(self,wp):
        pass
    
    
class MultNoise(Noise):
    '''
    ##################### PURPOSE: #############################
    Base-class for multiplicative noise; implements an efficient 
    evaluation of multiplicative noise with a realization of a
    Wiener process wp
    
    ############################################################
    '''
    def __init__(self,sigma): 
        self.sigma = sigma
        
    def eval_with_wp(self,v,wp):
        pass    
        

                              
        
             
                    
