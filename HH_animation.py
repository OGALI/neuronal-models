import numpy as np
import nice_colors as nc
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation


class HHAnimation(FuncAnimation):
    
    def __init__(self,spde_hh,
                 U_Jn,n_Jn,m_Jn,h_Jn,
                 plot_y_range): 
        
        # set SPDE
        self.spde_hh = spde_hh         
        
        # set grid
        self.grid = spde_hh.st_grid
        
        # resting potential
        self.rp = self.spde_hh.U.init_value.eq
        
        # set numerical solution for
        # potential
        self.U_Jn = U_Jn
        
        # ion channels, probabilities
        self.n_Jn = n_Jn
        self.m_Jn = m_Jn
        self.h_Jn = h_Jn
        
        # set wave speeds
        #self.wave_speed = wave_speed  
            
        # set plotting ranges
        self.plot_y_range = np.asarray(plot_y_range)
        
        # plotting range: x-spatial variable
        self.xlim = self.grid.spatial_grid.spatial_domain
        #self.tlim = self.grid.temporal_grid.temp_domain
        
        # plotting range: y-membrane potential
        self.yU_min = self.plot_y_range[1,0] 
        self.yU_max = self.plot_y_range[1,1] 
        
        # plotting range: y-ion channels
        self.y_ion_channels_min = self.plot_y_range[0,0]
        self.y_ion_channels_max = self.plot_y_range[0,1]
        
        # y-ranges
        self.U_range = (self.yU_min, self.yU_max)
        self.ion_channel_range = (self.y_ion_channels_min-0.2,
                                  self.y_ion_channels_max+1.)        
        
        # set figure
        self.fig = plt.figure()  
        
        # set upper sub-plot (opening probabilities ion channels)
        self.ax1 = self.fig.add_subplot(211, xlim = self.xlim,
                                        ylim = self.ion_channel_range)                     
        #self.ax1.set_xlabel('$x$')
        self.ax1.set_ylabel('$n(t),\,m(t),\,h(t)$')
        
        # initial data
        self.n_line, = self.ax1.plot([], [], 
                                     c=tuple(nc.colors[:,4]), 
                                     label=r'$n$')
        self.m_line, = self.ax1.plot([], [], 
                                     c=tuple(nc.colors[:,3]), 
                                     label=r'$m$')
        self.h_line, = self.ax1.plot([], [], 
                                     c=tuple(nc.colors[:,1]), 
                                     label=r'$h$')
        
        self.legend1 = self.ax1.legend(prop=dict(size=12))
        
        # set bottom sub-plot (membrane potential)
        self.ax2 = self.fig.add_subplot(212, xlim = self.xlim,
                                        ylim = self.U_range)
        self.ax2.axhline(self.rp, c='k', ls = ':', label=r"$resting\,pot.$")
        self.ax2.set_xlabel('$x,\,(mm)$')
        self.ax2.set_ylabel('$U(t,x),\,(mV)$')
        
        # initial data
        self.U_line, = self.ax2.plot([], [], c='r')
        
        # set title
        self.title = self.ax1.set_title("")
        
        # set legend
        self.legend2 = self.ax2.legend(prop=dict(size=12))
             
    # initiator; defines base frame for animation
    def init(self):
        self.title.set_text("")
        self.n_line.set_data([], [])
        self.m_line.set_data([], [])
        self.h_line.set_data([], [])
        self.U_line.set_data([], [])
        return (self.n_line, self.m_line,
                self.h_line, self.U_line, 
                self.title)

    # animator; updates animation sequentially
    # INPUT: n - frame number  
    def animate(self,n):     
        self.title.set_text("t = %.2f ms" % self.grid.temporal_grid.t[n])
        self.n_line.set_data(self.grid.spatial_grid.X_J, 
                             self.n_Jn[:,n])
        self.m_line.set_data(self.grid.spatial_grid.X_J, 
                             self.m_Jn[:,n])
        self.h_line.set_data(self.grid.spatial_grid.X_J, 
                             self.h_Jn[:,n])
        self.U_line.set_data(self.grid.spatial_grid.X_J, 
                             self.U_Jn[:,n])
        return (self.n_line, self.m_line, 
                self.h_line, self.U_line, 
                self.title)
    
    # run animation
    def run_save_animation(self):
        anim = animation.FuncAnimation(self.fig, 
                                       self.animate, 
                                       init_func = self.init, 
                                   frames = self.grid.temporal_grid.t.shape[0], 
                                   interval = 20, blit = True, repeat=False)
        
        anim.save('hh_det5.mp4', fps = None, 
                  extra_args=['-vcodec', 'libx264'])
        plt.show()
    
