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
    max_t_func = partial(max_epoch,max=13000)  # use partial functions and keyword args to create variable stop conditions (same function, multiple uses)

    t0=0
    y0 = [0,6378+500,0,np.sqrt(mu/(6378+500)),0,0]#spacecraft.y[-1]
    r0=y0[0:3]
    
    v0=y0[3:6]
    
    sc = []

    prop = Propagator()
    for i in range(50):
        sc.append(Spacecraft(t0,r0,v0))
        sc[i].propagate(prop, 10000, 1)
        sc[i].impulse_maneuver([-i/15,i/7.5,0])
        sc[i].propagate(prop,sc[i].t+19000,1)
        


    fig = sc[0].plot(show=False)

    for craft in sc[1:-1]:
        craft.add_plot(fig,show=False)
    
    sc[-1].add_plot(fig,show=True)

