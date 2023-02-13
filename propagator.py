import numpy as np
from scipy.integrate import ode


G = constants["g"]

class Propagator():
    

    def init(self, central_body, perturbations=None):
        
        self.central_body = central_body
        self.perturbations = perturbations
        self.solver = ode(self.EOM,) 
        
    def propagate(self, spacecraft, tf, dt):
        
        t0 = spacecraft.t[-1]

    def EOM(t,y):
        
        pos = y[:,0:3]
        vel = y[:,3:6]

        acc = -r * G / pow(r_mag,3)

        if j2 in lower(self.perturbations):

            acc += acc_j2
        
