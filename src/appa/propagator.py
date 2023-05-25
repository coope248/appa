import numpy as np
from scipy.integrate import ode

from .solarsystem import bodies
import appa.toolbox as tb

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
basic_propagator.frame = "ECLIPJ2000"

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
                "J2",
                "atmo_drag",
                ]
        if perturbation not in possible_perturbs:
            raise Exception('Perturbation type not valid')
            return
        if (perturbation == "J2") and ((self.frame!="J2000") or (self.central_body != "EARTH")):
            raise Exception("J2 only supported for Earth orbits using J2000 frame")
        if perturbation not in self.perturbations:
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
        t = np.ones((int(steps+1),1))*t0
        y = np.zeros((int(steps+1),6))
        t[0] = t0
        y[0] = y0
        self.solver.set_initial_value(y[0],t[0])
        i = 1

        stop = False
        if "low_thrust" in self.perturbations:
            self.thrust = spacecraft.thrust
            self.thrust_dim = np.size(self.thrust)
        if "atmo_drag" in self.perturbations:
            self.ballistic_coef = spacecraft.ballistic_coefficient
        while (self.solver.successful()) and (i <= steps) and (not stop):
            self.solver.integrate(self.solver.t+dt)
            t[i] = self.solver.t
            y[i] = self.solver.y
            i += 1
            if(bool(stop_cond)):
                stop = True in [f(self.solver.t,self.solver.y) for f in stop_cond]

        # if stopped, remove unused array indices
        if t0 in t[1:]:
            indices = np.where(t==t0)[0]
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
        mu = bodies[self.central_body]["mu"]
        cb_radius = bodies[self.central_body]["radius"]
        acc = -pos * mu/ pow(r_mag,3)
        if 'low_thrust' in self.perturbations:
            if  self.thrust_dim == 1:
                acc_lt = self.thrust * (vel / v_mag)

            elif self.thrust_dim == 3:
                #convert VNB reference frame to inertial
                v_hat = vel / np.linalg.norm(vel)
                h = np.cross(pos,vel)
                n_hat = h / np.linalg.norm(h)
                b_hat = np.cross(v_hat,n_hat)

                rot_mat = np.array([v_hat,n_hat,b_hat])
                acc_lt = np.dot(self.thrust,rot_mat)
            acc += acc_lt
        if "J2" in self.perturbations:
            j2 = bodies[self.central_body]["J2"]
            
            t = t = pos/r_mag*(5*pos[2]**2/r_mag**2)-np.dot(pos/r_mag,[[1,0,0],[0,1,0],[0,0,3]])

            acc_j2 = 1.5*j2*mu*cb_radius**2/r_mag**4*t
            acc += acc_j2

        if "atmo_drag" in self.perturbations:
            print(r_mag,r_mag-cb_radius)
            density = tb.get_atmo(self.central_body,r_mag - cb_radius)
            i_omega_cb = bodies[self.central_body]["rotation_rate_"+str(self.frame)]
            vel_rel = vel - np.cross(i_omega_cb,pos)
            v_r_mag = np.linalg.norm(vel_rel)
            e_hat_vr = vel_rel/v_r_mag
            acc_drag = -density/(2*self.ballistic_coef)*v_r_mag**2*e_hat_vr
            acc += acc_drag


    

        for body in self.bodies:
            body_mu = bodies[body]["mu"]
            r_cb_body = tb.ephemeris_getter(body,t,self.frame,self.central_body)[:3]
            r_sc_body = r_cb_body - pos
            r_sc_mag = np.linalg.norm(r_sc_body)
            r_cb_mag = np.linalg.norm(r_cb_body)
            
            acc_body = body_mu*(r_sc_body/r_sc_mag**3 - r_cb_body/r_cb_mag**3)
            acc += acc_body


        return [vel[0], vel[1], vel[2], acc[0], acc[1], acc[2]]
        
