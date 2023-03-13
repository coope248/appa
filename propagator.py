import numpy as np
from scipy.integrate import ode

from solarsystem import bodies

class Propagator():
    

    def __init__(self, central_body='Earth', perturbations=[]):
        
        self.central_body = central_body
        self.perturbations = perturbations
        self.solver = ode(self.EOM) 
        self.solver.set_integrator('dop853')
        
    def add_perturbation(self,perturbation):
        possible_perturbs = [
                'low_thrust',
                ]
        if perturbation not in possible_perturbs:
            raise Exception('Perturbation type not valid')
            return

        self.perturbations.append(perturbation)

        
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
        stop = False
        if 'low_thrust' in self.perturbations:
            self.thrust = spacecraft.thrust
        while (self.solver.successful()) and (i < steps) and (not stop):
            self.solver.integrate(self.solver.t+dt)
            t[i] = self.solver.t
            y[i] = self.solver.y
            i += 1
            if(bool(stop_cond)):
                stop = True in [f(self.solver.t,self.solver.y) for f in stop_cond]

        # if stopped, remove unused array indices
        if 0 in t[1:]:
            indices = np.where(t==0)[0]
            if indices[0]==0:
                first_zero = indices[1]
            else:
                first_zero = indices[0]
                
            t =  t[0:first_zero]
            y = y[~np.all(y==0,axis=1)]
            print(t.shape,y.shape)
        return t,y


    def EOM(self,t,y):
        
        pos = np.array(y[0:3])
        vel = np.array(y[3:6])
        
        r_mag = np.linalg.norm(pos)
        v_mag = np.linalg.norm(vel)
        acc = -pos * bodies[self.central_body]["mu"] / pow(r_mag,3)

        if 'low_thrust' in self.perturbations:
            acc_lt = self.thrust * (vel / v_mag)
            acc += acc_lt

        return [vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]]
        
