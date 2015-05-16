import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation


class NagumoAnimation(FuncAnimation):
    
    def __init__(self,spde_nagumo,v_Jn,wave_speed,plot_y_range): 
        
        # set Nagumo equation
        self.spde_nagumo = spde_nagumo         
        
        # set grid
        self.grid = spde_nagumo.st_grid
        
        # set numerical solution
        self.v_Jn = v_Jn
        
        # set wave speeds
        self.wave_speed = wave_speed  
            
        # set plotting ranges
        self.plot_y_range = np.asarray(plot_y_range)
        
        self.xlim = self.grid.spatial_grid.spatial_domain
        self.tlim = self.grid.temporal_grid.temp_domain
        
        self.yc_min = self.plot_y_range[0,0]
        self.yc_max = self.plot_y_range[0,1]
        
        self.yv_min = self.plot_y_range[1,0] 
        self.yv_max = self.plot_y_range[1,1]
        
        # y-ranges
        self.c_range = (self.yc_min -5. * (self.yc_max - self.yc_min),
                        self.yc_max + 20. * (self.yc_max - self.yc_min))        
        self.v_range = (self.yv_min, self.yv_max)
        
        # set figure
        self.fig = plt.figure()  
        
        # set upper sub-plot (for wave speed)
        self.ax1 = self.fig.add_subplot(211, xlim = self.tlim,
                                        ylim = self.c_range)             
        self.ax1.axhline(spde_nagumo.c_stoch, c='k', ls = ':',
                         label=r"$c=\sqrt{2b\nu}\left(\frac{1}{2} - a\right)$")
        self.ax1.set_xlabel('$t$')
        self.ax1.set_ylabel('$c(t)$')
        
        # initial data
        self.c_line, = self.ax1.plot([], [], c='r', 
                                     label='$c(t)$')
        
        # set bottom sub-plot (for TW)
        self.ax2 = self.fig.add_subplot(212, xlim = self.xlim,
                                        ylim = self.v_range)
        self.ax2.set_xlabel('$x$')
        self.ax2.set_ylabel('$v(t,x)$')
        
        # initial data
        self.wave_line, = self.ax2.plot([], [], c='r', 
                                        label='$v(t,x)$')
        
        # set title
        self.title = self.ax1.set_title("")
        
        # set legend
        self.legend = self.ax1.legend(prop=dict(size=12))
     
        
    # initiator; defines base frame for animation
    def init(self):
        self.title.set_text("")
        self.c_line.set_data([], [])
        self.wave_line.set_data([], [])
        return (self.c_line, 
                self.wave_line, 
                self.title)


    # animator; updates animation sequentially
    # INPUT: n - frame number  
    def animate(self,n):     
        self.title.set_text("t = %.2f" % self.grid.temporal_grid.t[n])
        self.c_line.set_data(self.grid.temporal_grid.t[:n], 
                             self.wave_speed[:n])
        self.wave_line.set_data(self.grid.spatial_grid.X_J, 
                                self.v_Jn[:,n])
        return (self.c_line, 
                self.wave_line, 
                self.title)
    
    
    # run animation
    def run_save_animation(self):
        anim = animation.FuncAnimation(self.fig, 
                                       self.animate, 
                                       init_func = self.init, 
                        frames = self.grid.temporal_grid.t.shape[0], 
                        interval = 20, blit = True, repeat=False)
        
        anim.save('wave_animation_nagumo6.mp4', 
                  fps = None, 
                  extra_args=['-vcodec', 'libx264'])
        plt.show()
    
