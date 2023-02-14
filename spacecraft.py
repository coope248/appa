import numpy as np

import plotly.express as px


class Spacecraft():

    def __init__(self, t0, r0, v0):
        self.t0 = t0
        self.y0 = y0
        self.ts = np.array([t0])
        self.ys = np.array([[y0,v0]])

    def propagate(self, propagator, tf, dt, stop_cond=None):
        propagator.propagate(self, tf, dt)

    def plot(self, show=True):
        # plot trajectory of spacecraft so far, 
        # returns plotly figure to allow multiple plots using add_plot method
        fig = px.line_3d(x = self.ys[:,0], 
                         y = self.ys[:,1], 
                         z = zself.ys[:,2], 
                        title = "orbit plot") 
        if show:
            fig.show
        return fig

    def add_plot(self, fig):
        #adds plot of spacecraft to previously created plotly figure
        ...
        

    def impulse_maneuver(self, delta_v):
        # add impulsive maneuver to trajectory, instantly changing the velocity vector
        
        #convert VNB reference frame to inertial
        ...
    
