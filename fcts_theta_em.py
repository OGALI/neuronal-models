import numpy as np
import scipy.sparse as sparse
from scipy.optimize import fsolve

import HH_drift_diff_coeff as hhd


def fiter(x,add_param,f):
    '''
    ################### PURPOSE : ###############
    composes function call for fsolve    
    
    ################### INPUT: ##################
    x               : variable
    add_param       : length-4 tuple
      add_param[0]  : theta
      add_param[1]  : delta_t
      add_param[2]  : right hand side
      add_param[3]  : tuple of additional 
                      parameters for f 
    f               : drift coefficient                 
    
    #############################################
    '''
    return (x-add_param[0]*add_param[1]*f(x,add_param[3])
            -add_param[2]) 

def theta_em_single_step(spde,est,add_param):
    f = spde.f.eval
    y = fsolve(fiter, est, (add_param,f))
    return y
    

#########################################################
############### complete implicit solver ###############
#########################################################
def compose_it_matrix(A,J,diff_param,delta_t,h,BC):
    # iteration matrix
    EE = sparse.eye(J) + diff_param * delta_t * A/h/h 
    return EE


class ImpEM(object):
    '''
    ###################### PURPOSE: ##########################
    Fully implicit scheme for discrete Hodgkin-Huxley system
    
    ##########################################################
    '''
    def __init__(self,hh_param,delta_t,h,J,BC,A):
        # PARAMETERS
        self.hh_param = hh_param    
                
        self.param_U = hh_param.param_U
        self.param_n = hh_param.param_n
        self.param_m = hh_param.param_m
        self.param_h = hh_param.param_h
        
        self.A = A
        self.J = J
        self.BC = BC
        
        self.delta_t = delta_t
        self.h = h
        
        self.EE = compose_it_matrix(self.A,
                                    self.J,
                                    self.param_U['diff_param'],
                                    self.delta_t,
                                    self.h,
                                    self.BC)
        
        # alpha_n, beta_n
        self.alpha_n = hhd.AlphaX(self.param_n)
        self.beta_n = hhd.BetaX(self.param_n)
        
        # alpha_m, beta_m
        self.alpha_m = hhd.AlphaX(self.param_m)
        self.beta_m = hhd.BetaX(self.param_m)
        
        # alpha_h, beta_h
        self.alpha_h = hhd.AlphaH(self.param_h)
        self.beta_h = hhd.BetaH(self.param_h)

    def g_n(self,U,n,m,h):
        return (n*(1.+(self.alpha_n.eval(U)+self.beta_n.eval(U))
                   *self.delta_t)-self.delta_t*self.alpha_n.eval(U))
            
    def g_m(self,U,n,m,h):
        return (m*(1.+(self.alpha_m.eval(U)+self.beta_m.eval(U))
                   *self.delta_t)-self.delta_t*self.alpha_m.eval(U))
    
    def g_h(self,U,n,m,h):
        return (h*(1.+(self.alpha_h.eval(U)+self.beta_h.eval(U))
                   *self.delta_t)-self.delta_t*self.alpha_h.eval(U))
    
    def f_U(self,U,n,m,h):
        return (self.param_U['gNa'] * (m**3) * h * (U-self.param_U['ENa'])
               + self.param_U['gK'] * (n**4) * (U-self.param_U['EK'])
               + self.param_U['gL'] * (U-self.param_U['EL']))
        
    def f(self,U,n,m,h):
        a = np.dot(self.EE.toarray(),U)
        b = self.delta_t * self.f_U(U,n,m,h)
        return a + b
                
    def f_iter(self,x,rs_U,rs_n,rs_m,rs_h):
        split = np.array_split(x,4)
        U = split[0]
        n = split[1]
        m = split[2]
        h = split[3]
        
        rs = np.hstack([rs_U,rs_n,rs_m,rs_h])
        
        U1 = self.f(U,n,m,h)
        n1 = self.g_n(U, n, m, h)
        m1 = self.g_m(U, n, m, h)
        h1 = self.g_h(U, n, m, h)
    
        ls = np.hstack([U1,n1,m1,h1])
        
        return rs - ls
    
    def imp_em_solve(self,init_guess,rs):
        solution = fsolve(self.f_iter,init_guess,rs)
        split = np.array_split(solution,4)
        U = split[0]
        n = split[1]
        m = split[2]
        h = split[3]
        return U,n,m,h 
    
    
    
