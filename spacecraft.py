import numpy as np

import plotly.express as px
import plotly.graph_objects as go


class Spacecraft():
    """
    A class that represents the physical aspects of a spacecraft in orbit. Holds trajectory, fuel, and plotting information

    ...

    Attributes
    ----------

    t : 
        current time value for craft
    y : 
        current state [x,y,z,vx,vy,vz]
    ts :
        array of all t values in spacecraft history
    ys : 
        array of all states corresponding to times in ts
    thrust :
        engine thrust of craft (used for finite burn calculations)


    Methods
    -------

    propagate(self, propagator, tf, dt, stop_cond=None):
        moves spacecraft forward in time according to EOMs defined by propagator 
    plot(show=True):
        plots the trajectory of the spacecraft, returns plotly plog figure
    add_plot(fig, show=True):
        adds the spacecraft's trajectory to the plotly figure passed as fig
    impulse_maneuver(delta_v):
        takes a delta_v in the VNB system and adds it to the current state, keeping t constant to simulate an instantaneous change in velocity of the craft

    """
    
    def __init__(self, t0, r0, v0):
        self.t = np.array(t0)
        self.y = np.array(r0+v0)
        self.ts = np.array([t0])
        self.ys = np.array([r0+v0])
        self.thrust = 0

    def propagate(self, propagator, tf, dt, stop_cond=None):
        t,y = propagator.propagate(self, tf, dt, stop_cond)
        self.ts = np.append(self.ts,t)
        self.ys = np.append(self.ys,y,0)
        self.t = t[-1]
        self.y = y[-1]


    def plot(self, show=True):
        # plot trajectory of spacecraft so far, 
        # returns plotly figure to allow multiple plots using add_plot method

        bound = np.absolute(self.ys).max() + 500
        xMax = [-bound,-bound,-bound,-bound,bound,bound,bound,bound]
        yMax = [-bound,-bound,bound,bound,-bound,-bound,bound,bound]
        zMax = [-bound,bound,-bound,bound,-bound,bound,-bound,bound]
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(x = xMax, 
                         y = yMax, 
                         z = zMax,
                        mode = 'markers',
                        marker=dict(
                            size=0.01,
                            opacity=0.01)))
        fig.add_trace(go.Scatter3d(
            x=self.ys[:,0],
            y=self.ys[:,1],
            z=self.ys[:,2],
                mode='lines',))

        if show:
            fig.show()
        return fig

    def add_plot(self, fig, show=True):
        #adds plot of spacecraft to previously created plotly figure
        fig.add_trace(go.Scatter3d(x = self.ys[:,0],
                                 y = self.ys[:,1],
                                 z = self.ys[:,2],
                                   mode = 'lines',
                                   line=dict(
                                            color='darkblue',
                                            width=2
                                            )   
                                   ))

        bound = np.absolute(self.ys).max() + 500
        xMax = [-bound,-bound,-bound,-bound,bound,bound,bound,bound]
        yMax = [-bound,-bound,bound,bound,-bound,-bound,bound,bound]
        zMax = [-bound,bound,-bound,bound,-bound,bound,-bound,bound]
        fig.update_traces(
                        x = xMax,
                        y = yMax,
                        z = zMax,
                        selector=dict(type="scatter3d", mode="markers"))
        if show:
            fig.show()
        

    def impulse_maneuver(self, delta_v):
        # add impulsive maneuver to trajectory, instantly changing the velocity vector
        delta_v = np.array(delta_v)
        r = self.y[0:3]
        v = self.y[3:6]

        #convert VNB reference frame to inertial
        v_hat = v / np.linalg.norm(v)
        h = np.cross(r,v)
        n_hat = h / np.linalg.norm(h)
        b_hat = np.cross(v_hat,n_hat)

        rot_mat = np.array([v_hat,n_hat,b_hat])
        delta_v = np.dot(delta_v,rot_mat)
        #add delta V to current velocity
        self.y += np.insert(delta_v,[0,0,0],0)
    
