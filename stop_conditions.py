import numpy as np
import toolbox as tb


def radius_stop(t, state, max_val=None, min_val=None, equals_val=None, accuracy=1):
    '''
    stop condition function to stop when orbital radius meets certain conditions

    Parameters:
    -----------

    t : float
        t value (required input for stop conditiion functions)
    state : array
        state of spacecraft, in for of [x, y, z, vx, vy, vz] (required input for stop conditions)
    max_val : float, optional
        maximum orbital radius value. If radius goes above this value, function returns true
    min_val : float, optional
        minimum orbital radius value. If radius goes below this value, function returns true
    equals_val : float, optional
        function returns true if radius equals this parameter (within accuracy)
    accuracy : float, optional
        accuracy for equals comparison (too small of value may cause equals comparison to never return true even if the radius crosses this value
    
    Returns:
    --------

    stop : bool
        returns true if propagator should stop, false otherwise

    Use:
    ----

    -Use functools partial to create a function that only takes t and state as input, setting max_val, min_val,
    and equals_val with partial

    ex: earth_collision = partial(radius_stop,min_val=earth_radius)
    prop
    '''
    gt_max = False
    lt_min = False
    eq_value = False
    if (max_val is None and min_val is None and equals_val is None):
       raise Exception("Must enter max, min, or equals value")
    
    r = np.linalg.norm(state[0:3])
    if max_val is not None:
        gt_max = r > max_val
    if min_val is not None:
        lt_min = r < min_val
    if equals_val is not None:
        eq_value = abs(r - equals_val) < accuracy
    stop = gt_max or lt_min or eq_value
    return stop


def kepler_stop(t, state, mu=398600.4,params=None, max_val=None, min_val=None, equals_val=None, accuracy=1):
    '''
    stop condition based on one or more keplerian orbital parameters

    Parameters:
    -----------
    
    t : float
        t value (required input for stop conditiion functions)
    state : array
        state of spacecraft, in for of [x, y, z, vx, vy, vz] (required input for stop conditions)
    mu : float, optional
        mu value for system, required for calculation of orbital parameters
    params : str
        string or list of strings corresponding to parameters to be checked
    max_val : float, optional
        value or list of maximum orbital parameter value. If any parameter goes above its corresponding value, function returns true
    min_val : float, optional
        value or list of minimum orbital parameter value. If parameter goes below its corresponding value, function returns true
    equals_val : float, optional
        value or list of function returns true if parameter equals this value (within accuracy)
    accuracy : float, optional
        value or list of values for accuracy for equals comparison (too small of value may cause equals comparison to never return true even if the parameter values cross this value
    '''

    kep_dict = tb.state2kep(state,mu)
    stop = False

    if bool(len(np.shape(params))):
        for i,param in enumerate(params):
            val = kep_dict[param]
            if max_val is not None:
                gt_max = val > max_val[i]
                stop |= gt_max
            if min_val is not None:
                lt_min = val < min_val[i]
                stop |= lt_min
            if equals_val is not None:
                eq_value = abs(val - equals_val[i]) < accuracy[i]
                stop |= eq_value
    else:
        val = kep_dict[params]
        if max_val is not None:
            gt_max = val > max_val
            stop |= gt_max
        if min_val is not None:
            lt_min = val < min_val
            stop |= lt_min
        if equals_val is not None:
            eq_value = abs(val - equals_val) < accuracy
            stop |= eq_value
        

    
    return stop            

