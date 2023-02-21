import numpy as np
import math


def kep2state(sma = 6800, ecc = 0, inc = 0, aop = 0, raan = 0, ta = 0):
    #calculate r and v at a given keplerian parameters
    ...

def modkep2state():
    #calculate r and v with given modified keplerian parameters
    ...

def state2kep(state):
    #calculate keplerian parameters from pos and vel (including C, r_p, r_a)
    x,y,z,vx,vy,vz = state




