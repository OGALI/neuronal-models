# taken from Ermentrout, Terman, 
# Mathematical Foundations of Neuroscience, 
# Chapter 1, p. 23f.

# PARAMETERS for membrane potential 
# UNITS: U [mV]; gX [mS/cm^3]; EX [mV]

phys_constants = {
                  'R':8.31451,
                  'F':96485.3,
                  }

param_K = {
           'C_in':400.,
           'C_out':20.,
           'z':1.           
           }

param_Na = {
            'C_in':50.,
            'C_out':440.,
            'z':1.
            }

param_Cl = {
            'C_in':40.,
            'C_out':560.,
            'z':-2.
            }

param_gen_env = {
                'T':293.2,
                'TC':20.,
                'T_base':6.3,
                'Q_10':3.,
                'param_K':param_K,
                'param_Na':param_Na,
                'param_Cl':param_Cl                 
                }

param_U = {
           'tau':1.,
           'diff_param':0.1,
           'gNa':120., 'gK':36., 'gL':0.3,
           'ENa':50., 'EK':-77., 'EL':-54.4
           }   
     
# PARAMETERS for opening probabilities
# Equation: n
param_n = {'type':'n',
           'ax_1':0.01, 'ax_2':0.1, 'A_x':55.,
           'bx_1':0.125, 'bx_2':1/80., 'B_x':65.
           }
    
# Equation: m
param_m = {'type':'m',
           'ax_1':0.1, 'ax_2':1/10., 'A_x':40.,
           'bx_1':4., 'bx_2':1/18., 'B_x':65.
           }     

# Equation: h    
param_h = {'type':'h',
           'ax_1':0.07, 'ax_2':1/20., 'A_x':65.,
           'bx_1':1., 'bx_2':1/10., 'B_x':35.
           }    


##############################################
############# HH-Parameter class #############
##############################################
class HHData:
    
    # PHYSICAL CONSTANTS
    phys_constants = phys_constants
    
    def __init__(self,
               param_gen_env,
               param_U,
               param_n,
               param_m,
               param_h,
               param_sigma,
               param_disc,
               BC):
                
        # ENVIRONMENTAL CONSTANTS
        self.param_gen_env = param_gen_env
        
        # MODEL parameters
        self.param_U = param_U
        self.param_n = param_n
        self.param_m = param_m
        self.param_h = param_h
        self.BC = BC 
        
        # NOISE parameters
        self.param_sigma = param_sigma
        
        # DISCRETIZATION parameters
        self.param_disc = param_disc

         
