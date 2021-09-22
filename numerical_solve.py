import matplotlib.pyplot as plt
import numpy as np

from sympy.solvers import solve
from sympy import Symbol

import time

def calc_L1(M1, M2, sep):
    # Solve the equation to get the location of the L1 point
    # print("Solving for the location of the L1 point...")
    
    r = Symbol('r', real=True)
    
    L1 = solve(
          (M2/(r*r))
        + (M1/(sep*sep))
        - ( r*(M1+M2) / (sep**3) )
        - M1/((sep-r)**2)
    )

    if len(L1) == 1:
        L1 = L1[0]
    else:
        print("Error! Multiple values of L1 found!")
        print(L1)

    # # L1 is currently measured from the secondary. Convert
    # # to be measured from the primary.
    # L1 = float(sep - L1)
    
    return L1

def f(M1, M2, sep):
    # Numerically get the location fo the L1 point
    r = np.linspace(0.9*sep, 0.1*sep, 1000)
    x = (M2/(r*r)) + (M1/(sep*sep)) - ( r*(M1+M2) / (sep**3) ) - M1/((sep-r)**2)

    # Distance from M2
    x = np.interp(0, x, r, left=np.nan, right=np.nan)

    return x

M1  = 2e30
M2  = 1e30
sep = 5.085e08

t1 = time.clock()
L1 = f(M1, M2, sep)
print("Numerically found L1 point:   {:2.3e} | {:2.3e}s".format(L1, time.clock()-t1))

t1 = time.clock()
L1 = calc_L1(M1, M2, sep)
print("Symbolically found L1 point:  {:2.3e} | {:2.3e}s".format(L1, time.clock()-t1))
