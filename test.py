import toolbox as tb
from propagator import Propagator
from spacecraft import Spacecraft
from functools import partial
import plotly.express as px

import numpy as np

from solarsystem import bodies

earth_radius = bodies["Earth"]["radius"]
earth_mu = bodies["Earth"]["mu"]
def max_epoch(t,state,max=20000):
    x,y,z,vx,vy,vz = state
    return t>max


if __name__ == "__main__":
    max_t_func = partial(max_epoch,max=13000)  # use partial functions and keyword args to create variable stop conditions (same function, multiple uses)

    t0=0
    y0 = [0,6378+500,0,np.sqrt(earth_mu/(earth_radius+500)),0,0]#spacecraft.y[-1]
    r0=y0[0:3]
    
    v0=y0[3:6]
    
    sc = []

    prop = Propagator()
    for i in range(5):
        sc.append(Spacecraft(t0,r0,v0))
        sc[i].propagate(prop, 10000, 10)
        prop.add_perturbation('low_thrust')
        sc[i].thrust = 0.00001*(i+1)
        #sc[i].impulse_maneuver([i/25,0,0])
        sc[i].propagate(prop,sc[i].t+11000,100)
        #sc[i].impulse_maneuver([1,0,0])
        sc[i].propagate(prop,sc[i].t+10000,100)

    for craft in sc:
    
        keps = [tb.state2kep(state) for state in craft.ys]
        es = [kep['ecc'] for kep in keps]

        eFig = px.line(x=craft.ts,y=es)
        eFig.show()
        


    fig = sc[0].plot(show=False)

    for craft in sc[1:-1]:
        craft.add_plot(fig,show=False)
    
    sc[-1].add_plot(fig,show=True)

