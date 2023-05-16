import numpy as np

import plotly.express as px
import plotly.graph_objects as go


class Spacecraft():
    """
    Class representing the physical aspects of a spacecraft in orbit. Holds trajectory, fuel, and plotting information

    ...

    Attributes
    ----------

    t : float 
        current time value for craft
    y : 3-vector
        current state [x,y,z,vx,vy,vz]
    ts : array
        array of all t values in spacecraft history
    ys : array of 3-vectors
        array of all states corresponding to times in ts
    thrust : float
        engine thrust of craft (used for finite burn calculations)
    color : tuple
        Tuple containing RGB values used for plot color
    name : str
        name of spacecraft mainly useful for plot legend


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
        '''
        Creates spacecraft object given initial state information

        Parameters:
        -----------

        t0 : float
            t0 (Epoch) of spacecraft given in seconds since J2000

        r0 : 3-vector
            position of spacecraft in km at t = t0 in J2000 Frame (not Ecliptic)

        v0 : 3-vector
            velocity of spacecraft in km/s at t = t0

        '''


        self.t = np.array(t0)
        self.y = np.array(r0+v0)
        self.ts = np.array([t0])
        self.ys = np.array([r0+v0])
        self.thrust = 0
        self.color = (0,0,0)
        self.name = None

    def propagate(self, propagator, tf, dt, stop_cond=None):
        '''
        propagates spacecraft trajcetory using given propagator
        
        Parameters:
        -----------

        propagator : orbit propagator object
            the propagator object used to calculate trajectory

        tf : float
            final time of trajectory. propagates until t = tf unless stop condition is met

        dt : float
            time between data points of state arrays (does not affect numeric integration dt which is taken care of by integrator)

        stop_cond : list, optional
            list of stop condition functions to pass to the propagator object. propagation will stop if any function in list returns true. functions should be in form of f(time, state) and return bool


        '''
        t,y = propagator.propagate(self, tf, dt, stop_cond)
        self.ts = np.append(self.ts,t)
        self.ys = np.append(self.ys,y,0)
        self.t = t[-1]
        self.y = y[-1]


    def plot(self, show=True, name=None, color=None):
        '''
        Plots all trajectory points in state arrays for spacecraft object

        Parameters:
        -----------

        show : bool, optional
            used to determine whether or not to show plot object after creating plotly figure object

        name : str, optional
            name used in legend, defaults to name of spacecraft if there is one

        color : tuple, optional
            tuple containing RGB color values to use for line plot

        Returns:
        -------

        fig : plotly figure object
            plotly figure that can be manipulated and/or passed to add_plot method (see plotly documentation for more information)

        '''
        if color == None:
            color = self.color
        if name == None:
            name = self.name
        bound = np.absolute(self.ys).max() + 500
        xMax = [-bound,-bound,-bound,-bound,bound,bound,bound,bound]
        yMax = [-bound,-bound,bound,bound,-bound,-bound,bound,bound]
        zMax = [-bound,bound,-bound,bound,-bound,bound,-bound,bound]
        fig = go.Figure()
        fig.add_trace(go.Scatter3d(
            x = xMax, 
            y = yMax, 
            z = zMax,
            showlegend=False,
            mode = 'markers',
            marker=dict(
                size=0.01,
                opacity=0.01)))
        fig.add_trace(go.Scatter3d(
            x=self.ys[:,0],
            y=self.ys[:,1],
            z=self.ys[:,2],
            customdata=self.ts,
            mode='lines',
            name = name,
            line=dict(
                color="rgb{}".format(color),
                width=2)
            hovertemplate='<br>x:%{x}<br>y:%{y}<br>z:%{z}<br>t:%{customdata}'))

        if show:
            fig.show()
        return fig

    def add_plot(self, fig, show=True, name=None, color=None):
        '''
        Adds trajectory of spacecraft to an existing plotly figure

        Parameters:
        -----------

        fig : plotly figure
            figure that trajectory plot is added to

        show : Bool, optional
            determines whether or not to show resulting figure after trajectory is added
        
        name : str, optional
            name used in legend, defaults to name of spacecraft if there is one

        color : tuple, optional
            tuple containing RGB color values to use for line plot

        '''
        
        if color == None:
            color = self.color
        if name == None:
            name = self.name
        bound_current = fig.data[0].x[-1]
        fig.add_trace(go.Scatter3d(
            x = self.ys[:,0],
            y = self.ys[:,1],
            z = self.ys[:,2],
            customdata = self.ts,
            mode = 'lines',
            name=name,
            line=dict(
                color="rgb{}".format(color),
                width=2)
            hovertemplate='<br>x:%{x}<br>y:%{y}<br>z:%{z}<br>t:%{customdata}'))

        bound = np.absolute(self.ys).max() + 500
        if bound_current > bound:
            bound = bound_current
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
        '''
        Adds an instantaneous delta-v to the spacecraft's state

        Parameters:
        -----------

        delta_v : 3-vector
            delta-V vector to add to spacecrafts current velocity (in VNB frame)
        '''

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
    
