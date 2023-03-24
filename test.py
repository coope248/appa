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

def hsv_to_rgb(h, s, v):
    h /= 360
    if s == 0.0: v*=255; return (v, v, v)
    i = int(h*6.) # XXX assume int() truncates!
    f = (h*6.)-i; p,q,t = int(255*(v*(1.-s))), int(255*(v*(1.-s*f))), int(255*(v*(1.-s*(1.-f)))); v*=255; i%=6
    if i == 0: return (v, t, p)
    if i == 1: return (q, v, p)
    if i == 2: return (p, v, t)
    if i == 3: return (p, q, v)
    if i == 4: return (t, p, v)
    if i == 5: return (v, p, q)

if __name__ == "__main__":
    tb.load_ephemeris()
    fig = tb.plot_body("MERCURY BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    tb.add_body_plot(fig,"VENUS BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    tb.add_body_plot(fig,"EARTH BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    tb.add_body_plot(fig,"MARS BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    tb.add_body_plot(fig,"JUPITER BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    tb.add_body_plot(fig,"PLUTO BARYCENTER",(0,100000000),10000,"ECLIPJ2000","SUN",False)
    fig.show()
    max_t_func = partial(max_epoch,max=13000)  # use partial functions and keyword args to create variable stop conditions (same function, multiple uses)
    
    moon0 = np.array(tb.ephemeris_getter("MOON",0,"J2000","EARTH"))
    t0=0
    y0 = [0,6378+450000,0,np.sqrt(earth_mu/(earth_radius+700000)),0,0]#spacecraft.y[-1]i
    y0 = [moon0[0],moon0[1]+10000,moon0[2],moon0[3]-0.7,moon0[4],moon0[5]]
    print(y0)
    r0=y0[0:3]
    v0=y0[3:6]
    
    sc = Spacecraft(t0,r0,v0)
    sc2 = Spacecraft(t0,r0,v0)
    prop = Propagator()
    prop.bodies = ['MOON']
    sc.propagate(prop,500000,100)
    prop.bodies = []
    sc2.propagate(prop,500000,100)

    fig2 = tb.plot_body("MOON",(0,2000000),1000,show=False)
    sc.add_plot(fig2, show=False)
    sc2.add_plot(fig2,show=True)

    scs = []
    y0 = [0,6378+500,0,np.sqrt(earth_mu/(earth_radius+500)),0,0]#spacecraft.y[-1]i
    r0=y0[0:3]
    v0=y0[3:6]
    N = 10
    for i in range(N):
        scs.append(Spacecraft(t0,r0,v0))
        scs[i].propagate(prop, 10000, 100)
        prop.add_perturbation('low_thrust')
        scs[i].thrust = 0.00001*(i+1)
        #scs[i].impulse_maneuver([i/25,0,0])
        scs[i].propagate(prop,scs[i].t+11000,100)
        scs[i].color = hsv_to_rgb(360 * i/N,1,1)

    #for craft in scs:
    
    #    keps = [tb.state2kep(state) for state in craft.ys]
    #    es = [kep['ecc'] for kep in keps]

    #    eFig = px.line(x=craft.ts,y=es)
    #    eFig.show()
        


    fig3 = scs[0].plot(show=False)

    for craft in scs[1:]:
        craft.add_plot(fig3,show=False)
    
    fig3.show()
