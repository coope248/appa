import numpy as np
import math


def kep2state(sma = 6800, ecc = 0, inc = 0, aop = 0, raan = 0, ta = 0, mu =398600.4):
    #calculate r and v at a given keplerian parameters
    # TODO: checks for hyperbolic orbit

    aol = ta+aop
    slr = sma * (1 - ecc**2)
    c3 = -mu / (2*sma)
    h = np.sqrt(mu * slr)

    r = slr / (1 + ecc * np.cos(ta))
    v = np.sqrt(2 * (c3 + mu/r))
    fpa = np.arccos( h / (r * v))

    if ta > np.pi:
        fpa *= -1

    v_rth = [v * np.sin(fpa), v * np.cos(fpa), 0]
    r_rth = [r, 0, 0]

    rot_mat = rth2eci(raan=raan, inc=inc, aol=aol)
    
    x,y,z = np.dot(rot_mat, r_rth)
    vx,vy,vz = np.dot(rot_mat, v_rth)

    return [x,y,z,vx,vy,vz]


def modkep2state(ra=6800, rp=6800, inc=0, aop=0, raan=0, ta=0, mu=398600.4):
    #calculate r and v with given modified keplerian parameters
    sma = (r_per+r_apo)/2
    ecc = (r_apo/sma) - 1
    return kep2state(sma, ecc, inc, aop, raan, ta, mu)

def state2kep(state, mu=398600.4):
    #calculate keplerian parameters from pos and vel (including C3, r_p, r_a, fpa)
    x,y,z,vx,vy,vz = state
    r = np.linalg.norm([x,y,z])
    r_hat = np.array([x,y,z])/r
    r_vec = r_hat*r
    v = np.linalg.norm([vx,vy,vz])
    v_hat = np.array([vx,vy,vz])/v
    v_vec = v_hat*v

    c3 = (v**2 / 2) - (mu/r)
    h = np.cross([x,y,z],[vx,vy,vz])
    h_mag = np.linalg.norm(h)
    h_hat = h/h_mag

    theta_hat = np.cross(h_hat,r_hat)

    slr = h_mag**2 / mu
    sma = -mu / (2*c3)
    ecc = np.sqrt(1 - (slr/sma))
    e_vec = np.cross(v_vec,h)/mu - r_hat
    e_hat = e_vec/ecc
    fpa = np.arcsin(np.dot(v_hat,r_hat)) 
    
    if np.dot(v_hat,r_hat) < 0:
        ta = 2*np.pi - np.arccos(np.dot(r_hat,e_hat))
    else:
        ta = np.arccos(np.dot(r_hat,e_hat))

    inc = np.arccos(h_hat[2])
    
    line_of_nodes = np.cross([0,0,1],h)
    n_mag = np.linalg.norm(line_of_nodes)
    n_hat = line_of_nodes/n_mag
    
    if line_of_nodes[1] < 0:
        raan = 2*np.pi - np.arccos(n_hat[0])
    else:
        raan = np.arccos(n_hat[0])

    if e_hat[2] < 0:
        aop = 2*np.pi - np.arccos(np.dot(n_hat,e_hat))
    else: 
        aop = np.arccos(np.dot(n_hat,e_hat))
    
    rp = sma*(1-ecc)
    ra = sma*(1+ecc)

    vp_mag = np.sqrt((1+ecc) * (mu/rp))
    
    va_mag = np.sqrt((1-ecc) * (mu/ra))
    
    dict = {'sma':sma,
            'ecc':ecc,
            'inc':inc,
            'aop':aop,
            'raan':raan,
            'ta':ta,
            'c3':c3,
            'fpa':fpa,
            'rp':rp,
            'ra':ra,
            'vp':vp_mag,
            'va':va_mag,
            }

    return dict
    




def rth2eci(raan=0, inc=0, aol=0):
    #calculate rotation matrix from r_hat, theta_hat, h_hat frame to ECI frame
    co,so = [np.cos(raan), np.sin(raan)]
    ci,si = [np.cos(inc), np.sin(inc)]
    ct,st = [np.cos(aol), np.sin(aol)]
    x = [co*ct - so*ci*st, -co*st - so*ci*ct, so*si]
    y = [so*ct + co*ci*st, -so*st + co*ci*ct, -co*si]
    z = [si*st, si*ct, ci]

    return np.array([x,y,z])

