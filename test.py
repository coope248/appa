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
sc.plot()
sc.impulse_maneuver([0,1,0])
sc.propagate(prop,sc.t+5000,1)
sc.plot()
print(np.absolute(sc.ys).max())
#t,y = prop.propagate(10000,1)
#x_vals = y[:,0]
#y_vals = y[:,1]
#z_vals = y[:,2]
#fig = px.line_3d(x = x_vals, 
#                         y = y_vals, 
#                         z = z_vals, 
#                        title = "orbit plot") 
#fig.show()

#sc.propagate(prop, 100, 1)
#delta_v = np.array([1,0,0])
#sc.impulsive_maneuver(delta_v)
