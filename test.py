import toolbox
from propagator import Propagator
from spacecraft import Spacecraft

import plotly.express as px

import numpy as np

mu = 398600.4

t0=0
y0 = [0,6378+500,0,np.sqrt(mu/(6378+500)),0,0]#spacecraft.y[-1]
r0=y0[0:3]
v0=y0[3:6]
sc = Spacecraft(t0,r0,v0)
prop = Propagator()
sc.propagate(prop, 10000, 1)
fig = sc.plot(show=False)
sc.impulse_maneuver([2,0,0])
sc.propagate(prop,sc.t+19000,1)
sc.add_plot(fig,show=True)

