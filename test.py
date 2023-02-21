import toolbox
from propagator import Propagator
from spacecraft import Spacecraft
from functools import partial
import plotly.express as px

import numpy as np

def max_epoch(t,state,max=20000):
    x,y,z,vx,vy,vz = state
    return t>max

if __name__ == "__main__":
    mu = 398600.4

    t0=0
    y0 = [0,6378+500,0,np.sqrt(mu/(6378+500)),0,0]#spacecraft.y[-1]
    r0=y0[0:3]
    v0=y0[3:6]
    sc = Spacecraft(t0,r0,v0)
    prop = Propagator()
    max_t_func = partial(max_epoch,max=13000)  # use partial functions and keyword args to create variable stop conditions (same function, multiple uses)
    sc.propagate(prop, 10000, 1, [max_t_func])
    fig = sc.plot(show=False)
    sc.impulse_maneuver([2,0,0])
    sc.propagate(prop,sc.t+19000,1,[max_t_func])
    sc.add_plot(fig,show=True)


