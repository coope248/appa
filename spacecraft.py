import numpy as np

import plotly.express as px


class Spacecraft():

    def __init__(self, t0, r0, v0):
        self.t = np.array(t0)
        self.y = np.array(r0+v0)
        self.ts = np.array([t0])
        self.ys = np.array([r0+v0])

    def propagate(self, propagator, tf, dt, stop_cond=None):
        t,y = propagator.propagate(self, tf, dt)
        self.ts = np.append(self.ts,t)
        self.ys = np.append(self.ys,y,0)
        self.t = t[-1]
        self.y = y[-1]


    def plot(self, show=True):
        # plot trajectory of spacecraft so far, 
        # returns plotly figure to allow multiple plots using add_plot method
        fig = px.line_3d(x = self.ys[:,0], 
                         y = self.ys[:,1], 
                         z = self.ys[:,2], 
                        title = "orbit plot") 
        if show:
            fig.show()
        return fig

    def add_plot(self, fig):
        #adds plot of spacecraft to previously created plotly figure
        ...


    def impulse_maneuver(self, delta_v):
        # add impulsive maneuver to trajectory, instantly changing the velocity vector
        delta_v = np.array([0,0,0]+delta_v)

        #convert VNB reference frame to inertial
        
        #add delta V to current velocity
        self.y += delta_v
    
