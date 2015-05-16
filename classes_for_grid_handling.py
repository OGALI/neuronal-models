import numpy as np

class SpatioTemporalGrid:
    def __init__(self,L,spatial_two_sided,J,T,N):
        
        # set spatial grid
        self.spatial_grid = SpatialGrid(L, J, spatial_two_sided)
        
        # set temporal grid
        self.temporal_grid = TemporalGrid(T, N)
        

class TemporalGrid:
    def __init__(self,T,N):
        self.set_temporal_grid(T, N)
        
    def set_temporal_grid(self,T,N):
        
        # array_like
        # length-(N+1) array of evenly spaced grid points        
        self.t = np.linspace(0,T,N+1)
        
        # time step
        self.delta_t = np.abs(self.t[1]-self.t[0])

        self.number_of_steps = N    
        self.temp_domain = (0,T)
     
        
class SpatialGrid:
    def __init__(self,L,J,spatial_two_sided):
        self.set_spatial_grid(L,J,spatial_two_sided)
        
    def set_spatial_grid(self,L,J,spatial_two_sided):
        self.spatial_two_sided = spatial_two_sided 
        if self.spatial_two_sided == True:    
            # array_like
            # length-(J+1) array of evenly spaced grid points
            self.X_J = np.linspace(-L,L,J+1)
            
            # tuple; start-, end point of domain
            self.spatial_domain = (-L,L)
        else:
            self.X_J = np.linspace(0,L,J+1)
            self.spatial_domain = (0,L)   
        
        # mesh size
        self.h = np.abs(self.X_J[1] - self.X_J[0]) 
        self.J = J
   
