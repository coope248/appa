import toolbox as tb
from propagator import Propagator
from spacecraft import Spacecraft
from functools import partial
import plotly.express as px

import numpy as np

from solarsystem import bodies

earth_radius = bodies["EARTH"]["radius"]
earth_mu = bodies["EARTH"]["mu"]
def max_epoch(t,state,max=20000):
    x,y,z,vx,vy,vz = state
    return t>max


if __name__ == "__main__":
    tb.load_ephemeris()
    max_t_func = partial(max_epoch,max=13000)  # use partial functions and keyword args to create variable stop conditions (same function, multiple uses)
    
    moon0 = np.array(tb.ephemeris_getter("MOON",0,"J2000","EARTH"))
    t0=0
    y0 = [0,6378+450000,0,np.sqrt(earth_mu/(earth_radius+700000)),0,0]#spacecraft.y[-1]i
    y0 = [moon0[0],moon0[1]+10000,moon0[2],moon0[3]+0.5,moon0[4],moon0[5]]
    r0=y0[0:3]
    v0=y0[3:6]
    
    sc = Spacecraft(t0,r0,v0)
    sc2 = Spacecraft(t0,r0,v0)
    prop = Propagator()
    sc.propagate(prop,5000000,100)
    prop.bodies = ['MOON']
    sc2.propagate(prop,5000000,100)

    fig = sc.plot(show=False)
    sc2.add_plot(fig,show=True)
    #for i in range(5):
    #    sc.append(Spacecraft(t0,r0,v0))
    #    sc[i].propagate(prop, 10000, 10)
    #    prop.add_perturbation('low_thrust')
    #    sc[i].thrust = 0.00001*(i+1)
    #    #sc[i].impulse_maneuver([i/25,0,0])
    #    sc[i].propagate(prop,sc[i].t+11000,100)
    #    #sc[i].impulse_maneuver([1,0,0])
    #    sc[i].propagate(prop,sc[i].t+10000,100)

    #for craft in sc:
    
    #    keps = [tb.state2kep(state) for state in craft.ys]
    #    es = [kep['ecc'] for kep in keps]

    #    eFig = px.line(x=craft.ts,y=es)
    #    eFig.show()
        


    #fig = sc[0].plot(show=False)

    #for craft in sc[1:-1]:
    #    craft.add_plot(fig,show=False)
    
    #sc[-1].add_plot(fig,show=True)

