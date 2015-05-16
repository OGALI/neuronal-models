import numpy as np

import classes_for_spde_handling as spde
    
        
class AddNoise2(spde.AddNoise):
    '''
    ##################### PURPOSE: ###########################
    Sub-class of AddNoise; implements an efficient 
    evaluation of additive noise with a Wiener process
    
    ##########################################################
    '''   
    def __init__(self,sigma,B):
        spde.AddNoise.__init__(self,sigma,B)
    
    # overwrite corresponding method in Base-class AddNoise
    def eval_with_wp(self,wp):
        return self.sigma*np.dot(self.B, wp)
    
    
class AddNoise1(spde.AddNoise):
    '''
    ##################### PURPOSE: ###########################
    Sub-class of AddNoise; as AddNoise0
    
    ##################### INPUT: #############################
    sigma   : parameter controlling noise amplitude
    J       : number of spatial grid points 
    
    ##########################################################
    '''
    def __init__(self,sigma,J):
        B = np.ones(J,)
        spde.AddNoise.__init__(self,sigma,B)
    
    def eval_with_wp(self,wp):
        return self.sigma * wp.sum() * self.B          
        
