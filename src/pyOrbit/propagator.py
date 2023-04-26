import numpy as np
from scipy.integrate import ode

from .solarsystem import bodies
import pyOrbit.toolbox as tb

class Propagator():
    '''
    Class that takes care of the actual integration/propagation of a trajectory for a spacecraft object
    
    ...

    Attributes:
    -----------

    central_body : str
        Name of body from solar system file that the propagator uses as the central body for gravitational acceleration

    perturbations : list
        list of perturbations to calculate and add to the central body's gravitaion in the EOM

    Methods:
    --------

    add_perturbation(perturbation):
        validates perturbation and adds it to the perturbation list attribute
    propagate(spacecraft, tf, dt, stop_cond=None):
        propagates passed spacecraft's trajectory
    EOM(t, y):
        System equations of motion, returns time derivative of state at given time and state
    '''
    

    def __init__(self, central_body='EARTH', perturbations=[]):
        '''
        Creates propagator object with given parameters

        Parameters:
        ----------
        
        central_body : str, optional
            The central body used for main gravitational acceleration in EOM
        perturbations : list, optional
            list of perturbations used in the EOM

        '''
        self.bodies = [] 
        self.frame = "J2000"
        self.central_body = central_body
        self.perturbations = perturbations
        self.solver = ode(self.EOM) 
        self.solver.set_integrator('dop853')
        
    def add_perturbation(self,perturbation):
        '''
        Validates given perturbation and adds to the propagator's internal list of perturbations

        Parameters:
        -----------

        perurbation : str
            name of perturbation to be validated and added to propagator

        Allowed Perturbations:
        ----------------------
        
        low_thrust: 
            adds a thrust acceleration (defined in spacecraft object) in direction of velocity vector

        '''
        possible_perturbs = [
                'low_thrust',
                ]
        if perturbation not in possible_perturbs:
            raise Exception('Perturbation type not valid')
            return

        self.perturbations.append(perturbation)

        
    def propagate(self, spacecraft, tf, dt, stop_cond=None):
        '''
        Performs the integration/propagation of the spacecraft's trajectory this function can also be called using the spacecraft's propagate object

        Parameters:
        -----------

        spacecraft : spacecraft object
            spacecraft to be propagated
        tf : float
            ending time value for propagation
        dt : float
            time difference between points in state outputs
        stop_cond : list, optional 
            list of stop condition functions. propagation will stop if any function in list returns true. functions should be in form of f(time, state) and return bool
            
        '''
        t0 = spacecraft.t
        y0 = spacecraft.y
        steps = np.ceil((tf-t0)/dt)
        t = np.zeros((int(steps+1),1))
        y = np.zeros((int(steps+1),6))
        t[0] = t0
        y[0] = y0
        
        self.solver.set_initial_value(y[0],t[0])
        i = 1
        stop = False
        if 'low_thrust' in self.perturbations:
            self.thrust = spacecraft.thrust
        while (self.solver.successful()) and (i <= steps) and (not stop):
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
        return t,y


    def EOM(self,t,y):
        '''
        Equations of motion for the system modeled by the propagator object

        Parameters:
        -----------

        t : float
            time of state
        y : 6-vector
            position in state space in form of [x, y, z, vx, vy, vz]

        Returns:
        --------

        y_prime : 6-vector
            dy/dt at given state and time essentially just [vx, vy, vz, ax, ay, az]

        '''
        pos = np.array(y[0:3])
        vel = np.array(y[3:6])
        
        r_mag = np.linalg.norm(pos)
        v_mag = np.linalg.norm(vel)
        acc = -pos * bodies[self.central_body]["mu"] / pow(r_mag,3)
        if 'low_thrust' in self.perturbations:
            acc_lt = self.thrust * (vel / v_mag)
            acc += acc_lt

        for body in self.bodies:
            r_body = tb.ephemeris_getter(body,t,self.frame,self.central_body)[0:3]
            r_body_sc = pos - r_body
            r_body_sc_mag = np.linalg.norm(r_body_sc)
            acc += -r_body_sc * (bodies[body]["mu"] / pow(r_body_sc_mag,3))


        return [vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]]
        
