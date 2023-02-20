import numpy as np
from scipy.integrate import ode


mu = 398602.545 #constants["g"]

class Propagator():
    

    def __init__(self, central_body='earth', perturbations=[]):
        
        self.central_body = central_body
        self.perturbations = perturbations
        self.solver = ode(self.EOM) 
        self.solver.set_integrator('dop853')
        

        
    def propagate(self, spacecraft, tf, dt, stop_cond=None):
        
        t0 = spacecraft.t
        y0 = spacecraft.y
        steps = np.ceil((tf-t0)/dt)
        t = np.zeros((int(steps),1))
        y = np.zeros((int(steps),6))
        t[0] = t0
        y[0] = y0
        
        self.solver.set_initial_value(y[0],t[0])
        
        i = 1
        while (self.solver.successful()) and (i < steps) and (not stop):
            self.solver.integrate(self.solver.t+dt)
            t[i] = self.solver.t
            y[i] = self.solver.y
            i += 1
            if(bool(stop_cond)):
                stop = not False in [bool(f(solver.t,solver.y)) for f in stop_cond]

        
        return t,y
    def EOM(self,t,y):
        
        pos = np.array(y[0:3])
        vel = np.array(y[3:6])
        
        r_mag = np.linalg.norm(pos)

        acc = -pos * mu / pow(r_mag,3)

        if 'j2' in [p.lower() for p in self.perturbations]:
            acc += acc_j2

        return [vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]]
        
