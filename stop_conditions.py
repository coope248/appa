import numpy as np


def radius_stop(t,state,max_val=None,min_val=None,equals_val=None, accuracy=1):
    '''

    '''

    if (max_val is None and min_val is None and equals_val is None):
       raise Exception("Must enter max, min, or equals value")
    
    r = np.linalg.norm(state[0:3])
    if max_val is not None:
        gt_max = r > max_val
        print("gt:",gt_max)
    if min_val is not None:
        lt_min = r < min_val
        print("lt:",lt_min)
    if equals_val is not None:
        eq_value = abs(r - equals_val) < accuracy
        print("eq:",eq_value)


