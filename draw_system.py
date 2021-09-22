import numpy as np
from scipy.constants import h, c, k, G

from sympy.solvers import solve
from sympy import Symbol

import matplotlib.pyplot as plt
import matplotlib.patches as patch

from colour_system import cs_hdtv
cs = cs_hdtv


def mask(mask, X, Y, target):
    # Get info
    x0 = mask['x0']
    y0 = mask['y0']
    radius = mask['radius']
    
    # Radius from object at each point
    r = np.sqrt((X - x0)**2 + (Y - y0)**2)
    # Create mask
    outside = r > radius

    # Deploy mask
    target *= outside

    return target

def planck(lam, T):
    """ Returns the spectral radiance of a black body at temperature T.

    Returns the spectral radiance, B(lam, T), in W.sr-1.m-2 of a black body
    at temperature T (in K) at a wavelength lam (in nm), using Planck's law.

    """

    lam_m = lam / 1.e9
    fac = h*c/lam_m/k/T
    B = 2*h*c**2/lam_m**5 / (np.exp(fac) - 1)
    return B

def separation(M1, M2, period):
    # Compute separation from period and mass, assuming circular orbit
    separation = G * (M1 + M2) * (period**2)
    separation = separation / (4 * np.pi * np.pi)
    separation = separation ** (1/3) # in metres
    return separation

def roche(M1, M2, period, x, y, z=0.0):
    # Mass
    mu = M2/(M1+M2)
    # print("M1: {:.3e} --- M2: {:.3e}".format(M1, M2))
    
    a = separation(M1, M2, period)
    # print("With a period of {:.3e} and a separation of {:.3e} m".format(period, a))

    # Point distance from M1
    r1 = np.linalg.norm(
        np.array([x, y, z])
        )
    # Point distance from M2
    r2 = np.linalg.norm(
        np.array([(x-a), y, z])
        )

    # Angular velocity
    w = 2*np.pi / period #rad/s
    # print("Angular velocity: {:.3e} rad/s".format(w))

    # Get the angular momentum contribution
    V = (x-(mu*a))**2 + (y**2)
    V *= 0.5 * w**2
    # Add in gravitational potential
    V = V + (G*M1/r1) + (G*M2/r2)

    # Potential is negative
    V = - V

    return V

def calc_r_roche(M1, M2, separation):
    # Calculate the roce radius of the secondary (it fills this volume)
    q = M1/M2

    r = 0.49 * (q**(2/3))
    r /= (0.6*(q**(2/3))) + np.log(1 + (q**(1/3)))

    r *= separation

    return r

def calc_r_circ(R_L1, M1, period):
    # Calculate the circularisation radius of material ejected from the L1 point
    r = 4*np.pi*np.pi * (R_L1**4)
    r /= G * M1 * period * period

    return r

def OLD_calc_L1(M1, M2, sep):
    # Solve the equation to get the location of the L1 point
    print("Solving for the location of the L1 point...")
    
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

    # L1 is currently measured from the secondary. Convert
    # to be measured from the primary.
    L1 = float(sep - L1)
    
    print("Done! L1 is {:.2e} m from M1".format(L1))

    return L1

def calc_L1(M1, M2, sep):
    # Numerically get the location fo the L1 point
    r = np.linspace(0.9*sep, 0.1*sep, 1000)
    x = (M2/(r*r)) + (M1/(sep*sep)) - ( r*(M1+M2) / (sep**3) ) - M1/((sep-r)**2)

    # Distance from M2
    x = np.interp(0, x, r, left=np.nan, right=np.nan)

    return sep - x

def calc_r_WD(M):
    R_sun = 6.957e8
    M_sun = 1.98e30

    r = (M_sun/M)**(1/3)
    r *= R_sun * 0.010

    return r

def calc_E_lambda(l, T):
    E = 8*np.pi*h*c
    E /= (l**5) * (np.exp(h*c/(l*k*T)) - 1)

    return E

################# System properties #################
P        = 2   *60*60         # s
M1       = 1e30               # kg
M2       = 0.5e30             # kg
r_disc   = 0.2                # Scaled to separation
T_WD     = 30000              # Kelvin
T_sec    = 3000               # Kelvin
T_disc   = 2000               # Kelvin
T_BS     = 7000               # Kelvin
BS_az    = np.deg2rad(45)     # bright spot azimuth
#####################################################


# Calculate extra parameters
sep = separation(M1, M2, P)
r_disc *= sep
r_wd = calc_r_WD(M1)
r_sec = calc_r_roche(M1, M2, sep)
L1 = calc_L1(M1, M2, sep)
r_circ = calc_r_circ(L1, M1, P)


### Calculate colours of the bodies from temperatures
# The grid of visible wavelengths
lam = np.arange(380., 781., 5)

# Calculate the black body spectrum and the HTML hex RGB colour string
#WD
spec_WD = planck(lam, T_WD)
html_rgb_WD = cs.spec_to_rgb(spec_WD, out_fmt='html')
#Secondary
spec_sec = planck(lam, T_sec)
html_rgb_sec = cs.spec_to_rgb(spec_sec, out_fmt='html')
#Disc
spec_disc = planck(lam, T_disc)
html_rgb_disc = cs.spec_to_rgb(spec_disc, out_fmt='html')
#BS
spec_BS = planck(lam, T_BS)
html_rgb_BS = cs.spec_to_rgb(spec_BS, out_fmt='html')


# Report the physical parameters we calculated
print("Primary radius:           {:.3e} m".format(r_wd))
print("Secondary radius:         {:.3e} m".format(r_sec))
print("Separation:               {:.3e} m".format(sep))
print("L1:                       {:.3e} m from primary".format(L1))
print("Circularisation radius:   {:.3e} m".format(r_circ))
print("Disc radius:              {:.3e} m".format(r_disc))


# Initialise the plotting area
fig, ax = plt.subplots()
ax.axis('off')
ax.set_aspect('equal')
grey = (0.15, 0.15, 0.15)
ax.set_facecolor(grey)
fig.patch.set_facecolor(grey)

# Create the circle filled by the primary at coordinates 0,0
patch_WD = patch.Circle((0,0), r_wd, color=html_rgb_WD, zorder=0.)
ax.add_patch(patch_WD)

# Create the circle filled by the secondary
patch_sec = patch.Circle((sep, 0), r_sec, color=html_rgb_sec, zorder=0.)
ax.add_patch(patch_sec)
ax.scatter(sep, 0, color='grey', marker='x', s=7, zorder=10.)

# Create the annulus filled by the disc. Extends from R_disc in to R_inner_disc, or to 0?
# Placeholder value of 5*r_wd
r_disc_inner = 5*r_wd

n, radii = 100, [r_disc_inner, r_disc]
theta = np.linspace(0, 2*np.pi, n, endpoint=True)
xs = np.outer(radii, np.cos(theta))
ys = np.outer(radii, np.sin(theta))

# in order to have a closed area, the circles
# should be traversed in opposite directions
xs[1,:] = xs[1,::-1]
ys[1,:] = ys[1,::-1]

ax.fill(np.ravel(xs), np.ravel(ys), facecolor=html_rgb_disc, edgecolor=None, zorder=5.)


# Plot a dashed circle at r_circ
patch_r_circ = patch.Circle((0,0), r_circ, zorder=10., fill=False, edgecolor='black', linestyle='--')
ax.add_patch(patch_r_circ)

# Put a cross at the L1 point
ax.scatter(L1, 0, color='grey', marker='x', zorder=10.)

# Plot the bright spot at the point where a line intersecting the x axis is 
# tangent to the disc outer edge, at an angle BS_az
theta = BS_az
BS_x = r_disc * np.cos(theta)
BS_y = r_disc * np.sin(theta)
plt.scatter(BS_x, BS_y, s=20, color=html_rgb_BS, zorder=10.)

# Plot the BS strip as a colourbar, showing the light output as a function of space

# Save and show the figure
plt.savefig('system.svg')
plt.show()