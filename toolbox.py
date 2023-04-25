import numpy as np
import math
import os

import plotly.graph_objects as go
import plotly.express as px
import spiceypy as spice
from spice_data.create_mk import create_mk

SPICE_PATH = os.path.join(os.path.dirname(__file__),"spice_data")

def kep2state(sma = 6800, ecc = 0, inc = 0, aop = 0, raan = 0, ta = 0, mu =398600.4):
    '''
    Function to calculate state vector (r and v) at given keplerian parameters

    Parameters:
    -----------

    sma : float, optional
        semi major axis of keplerian orbit
    ecc : float, optional
        eccentricity of keplerian orbit
    inc : float, optional
        inclination of orbit
    aop : float, optional
        argument of periapsis
    raan : float, optional
        right ascension of the ascending node
    ta : float, optional
        true anomoly of object in radians
    mu : float, optional
        gravitational parameter of central body for keplerian orbit

    Returns:
    --------

    state : array
        state vector in form of [x, y, z, vx, vy, vz] corresponding to same orbit and position that is input

    '''

    #TODO: checks for hyperbolic orbit
    aol = ta+aop
    slr = sma * (1 - ecc**2)
    c3 = -mu / (2*sma)
    h = math.sqrt(mu * slr)

    r = slr / (1 + ecc * math.cos(ta))
    v = math.sqrt(2 * (c3 + mu/r))
    fpa = np.arccos( max(-1,min(1,h / (r * v))))

    if ta > np.pi:
        fpa *= -1

    v_rth = [v * math.sin(fpa), v * math.cos(fpa), 0]
    r_rth = [r, 0, 0]

    rot_mat = rth2eci(raan=raan, inc=inc, aol=aol)
    
    x,y,z = np.dot(rot_mat, r_rth)
    vx,vy,vz = np.dot(rot_mat, v_rth)

    return [x,y,z,vx,vy,vz]


def modkep2state(ra=6800, rp=6800, inc=0, aop=0, raan=0, ta=0, mu=398600.4):
    '''
    Function to calculate state vector (r and v) at given modified keplerian parameters (based on GMAT modified keplerian option)

    Parameters:
    -----------

    ra : float, optional
        distance at apoapsis of orbit
    rp : float, optional
        distance at periapsis of orbit
    inc : float, optional
        inclination of orbit
    aop : float, optional
        argument of periapsis
    raan : float, optional
        right ascension of the ascending node
    ta : float, optional
        true anomoly of object in radians
    mu : float, optional
        gravitational parameter of central body for keplerian orbit

    Returns:
    --------

    state : array
        state vector in form of [x, y, z, vx, vy, vz] corresponding to same orbit and position that is input
    '''
    sma = (r_per+r_apo)/2
    ecc = (r_apo/sma) - 1
    return kep2state(sma, ecc, inc, aop, raan, ta, mu)

def state2kep(state, mu=398600.4):
    '''
    calculates orbital parameters given the current state vector

    Parameters: 
    -----------

    state : array
        state vector in form of [x, y, z, vx, vy, vz] corresponding to same orbit and position that is input
    mu : float, optional 
        gravitational parameter of central body for keplerian orbit

    Returns:
    --------

    param_dict : dictionary
        dictionary of orbital parameters
        
    Available orbital parameters:
    -----------------------------
    
    ra : float
        distance at apoapsis of orbit
    rp : float
        distance at periapsis of orbit
    inc : float
        inclination of orbit
    aop : float
        argument of periapsis
    raan : float
        right ascension of the ascending node
    ta : float
        true anomoly of object in radians
    c3 : float
    fpa : float
    rp : float
    ra : float
    vp : float
    va : float
        

    '''
    x,y,z,vx,vy,vz = state
    r = np.linalg.norm([x,y,z])
    r_hat = np.array([x,y,z])/r
    r_vec = r_hat*r
    v = np.linalg.norm([vx,vy,vz])
    v_hat = np.array([vx,vy,vz])/v
    v_vec = v_hat*v

    c3 = (v**2 / 2) - (mu/r)
    h = np.cross([x,y,z],[vx,vy,vz])
    h_mag = np.linalg.norm(h)
    h_hat = h/h_mag

    theta_hat = np.cross(h_hat,r_hat)

    slr = h_mag**2 / mu
    sma = -mu / (2*c3)
    ecc = np.sqrt(1 - (slr/sma))
    e_vec = np.cross(v_vec,h)/mu - r_hat
    if ecc != 0:
        e_hat = e_vec/ecc
    else:
        e_hat = np.array([0,0,0])
    fpa = np.arcsin(np.dot(v_hat,r_hat)) 
    
    if np.dot(v_hat,r_hat) < 0:
        ta = 2*np.pi - np.arccos(np.dot(r_hat,e_hat))
    else:
        ta = np.arccos(max(-1,min(1,np.dot(r_hat,e_hat))))

    inc = np.arccos(h_hat[2])
    
    line_of_nodes = np.cross([0,0,1],h)
    n_mag = np.linalg.norm(line_of_nodes)
    if n_mag != 0:
        n_hat = line_of_nodes/n_mag
    else:
        n_hat = np.array([0,0,0])
    if line_of_nodes[1] < 0:
        raan = 2*np.pi - np.arccos(n_hat[0])
    else:
        raan = np.arccos(n_hat[0])

    if e_hat[2] < 0:
        aop = 2*np.pi - np.arccos(np.dot(n_hat,e_hat))
    else: 
        aop = np.arccos(np.dot(n_hat,e_hat))
    
    rp = sma*(1-ecc)
    ra = sma*(1+ecc)

    vp_mag = np.sqrt((1+ecc) * (mu/rp))
    
    va_mag = np.sqrt((1-ecc) * (mu/ra))
    
    param_dict = {'sma':sma,
            'ecc':ecc,
            'inc':inc,
            'aop':aop,
            'raan':raan,
            'ta':ta,
            'c3':c3,
            'fpa':fpa,
            'rp':rp,
            'ra':ra,
            'vp':vp_mag,
            'va':va_mag,
            }

    return param_dict
    




def rth2eci(raan=0, inc=0, aol=0):
    #calculate rotation matrix from r_hat, theta_hat, h_hat frame to ECI frame
    co,so = [math.cos(raan), math.sin(raan)]
    ci,si = [math.cos(inc), math.sin(inc)]
    ct,st = [math.cos(aol), math.sin(aol)]
    x = [co*ct - so*ci*st, -co*st - so*ci*ct, so*si]
    y = [so*ct + co*ci*st, -so*st + co*ci*ct, -co*si]
    z = [si*st, si*ct, ci]

    return np.array([x,y,z])


def spice_object_getter(file, printopt=False):
    objects = spice.spkobj(file)
    ids,names,tcs_s,tcs_pretty = [],[],[],[]

    if printopt:
        print("Objects retrived from {}".format(file))

    n=0 
    for obj in objects:
        ids.append(obj)
    
        tc_s = spice.wnfetd(spice.spkcov(file,ids[n]),n)
        tcs_s.append(tc_s)
        tcs_pretty.append([spice.timout(f,"YYYY MON DD HR:MM:SC.### (TDB) ::TDB" ) for f in tc_s])

        try:
            name = spice.bodc2n(obj)
        except:
            name = "UNKOWN"

        names.append(name)

        if printopt:
            print("ID: {0:3}\t\tName: {1:35}\tTime Cov: {2} --> {3}".format(ids[-1],names[-1],tcs_pretty[-1][0],tcs_pretty[-1][1]))

        

    return ids,names,tcs_s,tcs_pretty

def ephemeris_getter(target,times,frame,observer):
    return np.array(spice.spkezr(target,times,frame,'NONE',observer)[0])

def tc_array(tcs,n_steps):
    arr = np.zeros((n_steps,1))
    arr[:,0] = np.linspace(tcs[0],tcs[1],n_steps)
    return arr
    
def load_ephemeris():
    create_mk()
    spice.furnsh(os.path.join(SPICE_PATH,"ss_kernel.mk"))
    ids,names,tcs_s,tcs_pr = spice_object_getter(os.path.join(SPICE_PATH,"de432s.bsp"),True)
    return ids,names,tcs_s,tcs_pr


def plot_body(body,t, steps, frame="J2000",observer="EARTH", show=True, color=None):
    '''
    Plots all trajectory points in state arrays for spacecraft object

    Parameters:
    -----------
    body : str 
        body to plot
    t : tuple
        timespan to plot body (t0,tf)

    frame : str, optional
        frame to plot bodies trajectory defaults to J2000 frame

    observer : str, optional
        observer of trajectory (center of reference) defaults to EARTH

    show : bool, optional
        used to determine whether or not to show plot object after creating plotly figure object

    Returns:
    -------

    fig : plotly figure object
        plotly figure that can be manipulated and/or passed to add_plot method (see plotly documentation for more information)

    bound : float
        bounds of resulting plot figure
    '''
    
    if color == None:
        color = (0,0,0)

    ys = ephemeris_getter(body,tc_array(t,steps),frame,observer)
    bound = np.absolute(ys).max() + 500
    xMax = [-bound,-bound,-bound,-bound,bound,bound,bound,bound]
    yMax = [-bound,-bound,bound,bound,-bound,-bound,bound,bound]
    zMax = [-bound,bound,-bound,bound,-bound,bound,-bound,bound]
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(x = xMax, 
                               y = yMax, 
                               z = zMax,
                               mode = 'markers',
                               showlegend=False,
                               marker=dict(
                                   size=0.01,
                                   opacity=0.01)))
    fig.add_trace(go.Scatter3d(x=ys[:,0],
                               y=ys[:,1],
                               z=ys[:,2],
                               mode='lines',
                               line=dict(color="rgb{}".format(color),
                                         width=2)))

    if show:
        fig.show()
    return fig
    
def add_body_plot(fig,body, t, steps, frame="J2000",observer="EARTH",show=True,color=None):

    '''
    Adds trajectory of spacecraft to an existing plotly figure
        
    Parameters:
    -----------

    fig : plotly figure
        figure that trajectory plot is added to

    body : str
        body to add to figure

    t : tuple
        timespan to plot body (t0,tf)

    frame : str, optional
        frame to plot bodies trajectory defaults to J2000 frame

    observer : str, optional
        observer of trajectory (center of reference) defaults to EARTH

    show : Bool, optional
        determines whether or not to show resulting figure after trajectory is added

    '''
    
    if color == None:
        color = (0,0,0)

    ys = ephemeris_getter(body,tc_array(t,steps),frame,observer)
    bound_current = fig.data[0].x[-1]
    fig.add_trace(go.Scatter3d(x = ys[:,0],
                               y = ys[:,1],
                               z = ys[:,2],
                               mode = 'lines',
                               line=dict(color="rgb{}".format(color),
                                         width=2)))
        
    bound = np.absolute(ys).max() + 500
    if bound_current > bound:
        bound = bound_current
    xMax = [-bound,-bound,-bound,-bound,bound,bound,bound,bound]
    yMax = [-bound,-bound,bound,bound,-bound,-bound,bound,bound]
    zMax = [-bound,bound,-bound,bound,-bound,bound,-bound,bound]
    fig.update_traces(x = xMax,
                      y = yMax,
                      z = zMax,
                      selector=dict(type="scatter3d", mode="markers"))
    if show:
        fig.show()
        




