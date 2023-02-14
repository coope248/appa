import toolbox
from propagator import Propagator
from spacecraft import Spacecraft

import plotly.express as px

import numpy as np

#sc = Spacecraft()
prop = Propagator()
t,y = prop.propagate(10000,1)
x_vals = y[:,0]
y_vals = y[:,1]
z_vals = y[:,2]
fig = px.line_3d(x = x_vals, 
                         y = y_vals, 
                         z = z_vals, 
                        title = "orbit plot") 
fig.show()

print(len(y))
#sc.propagate(prop, 100, 1)
#delta_v = np.array([1,0,0])
#sc.impulsive_maneuver(delta_v)
