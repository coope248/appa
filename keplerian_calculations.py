import numpy as np
import math


def rv_calc(SMA = 6800, ECC = 0, INC = 0, AOP = 0, RAAN = 0, TA = 0):
    #calculate r and v at a given true anomaly
    ...

def hohmann_transfer(*args):
    nargs = len(args)
    assert((nargs % 2 == 0) && (nargs < 6)) # function requires 2-6 arguments describing parameters for each orbit
    if(nargs == 2):
        sma1 = args[0]
        sma2 = args[1]
        ecc1 = 0
        ecc2 = 0
        aop1 = 0
        aop2 = 0
    elif(nargs == 4):
        sma1 = args[0]
        ecc1 = args[1]
        sma2 = args[2]
        ecc2 = args[3]
        aop1 = 0
        aop2 = 0
    else:
        sma1 = args[0]
        ecc1 = args[1]
        aop1 = args[2]
        sma2 = args[3]
        ecc2 = args[4]
        aop2 = args[5]

    


	
	
	
def main():
    hohmann_transfer(2)
  
if __name__ == "__main__":
	main()  
 
