import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patch


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

def separation(M1, M2, period):
    G = 6.67e-11 # kg m/s**2
    # Compute separation from period and mass, assuming circular orbit
    separation = G * (M1 + M2) * (period**2)
    separation = separation / (4 * np.pi * np.pi)
    separation = separation ** (1/3) # in metres
    return separation

def roche(M1, M2, period, x, y, z=0.0):
    G = 6.67e-11 # kg m/s**2
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

M1  = 1e30         #kg
M2  = 0.5e30       #kg
R1 = 0.2           # Solar radii
Rdisc = 0.2 # Scaled to separation
BS_az = np.deg2rad(45)
P   = 2.5 *60*60   #s
sep = separation(M1, M2, P)
Rdisc *= sep


mask1 = {
    'x0': 0.0,
    'y0': 0.0,
    'radius': 0.2*sep # 6.37e7 *2.3
}
mask2 = {
    'x0': sep,
    'y0': 0.0,
    'radius': 0.2*sep # 0.3*6.37e7 *2
}

# Make a meshgrid
x = np.linspace(-2*sep, 2.5*sep, 100)
y = np.linspace(-2*sep, 2*sep, 100)

X, Y = np.meshgrid(x, y)

# Compute roche potential
V = roche(M1, M2, P, X, Y)

patches = []
for m in [mask1, mask2]:
    x0 = m['x0']
    y0 = m['y0']
    radius = m['radius']
    patches.append(patch.Circle((x0,y0), radius, fill=True, facecolor='blue'))


# from sympy.solvers import solve
# from sympy import Symbol

# print("Solving for the location of the L1 point...")
# r = Symbol('r', real=True)
# # L1 = solve(
# #     (M2/(r*r))
# #     + (M1/(sep*sep))
# #     - ( r * (M1+M2) / (sep**3) )
# #     - M1/((sep-r)**2)
# # )
# L1 = [3.37e+08]
# L1 = L1[0]
# L1 = float(sep-L1)
# print("Done! L1 is {:.2e} m from M1".format(L1))

# L1 = roche(M1, M2, P, L1, 0.0, 0.0)
# Contours = np.linspace(L1, np.nanmax(V), 15)
# Contours.sort()

# import scipy.interpolate as interpolate
# lobe = interpolate.griddata((x, y), V, L1)
# print(lobe)


#### PLOTTING
# Mask out near the bodies
V = mask(mask1, X, Y, V)
V = mask(mask2, X, Y, V)
V[V==0] = np.nan # ragged edges!



from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_wireframe(X, Y, V, rcount=30, ccount=30)


for i, row in enumerate(V):
    for j in range(1, len(row)-1):
        if ((row[j-1] < row[j] and row[j+1] < row[j]) or 
            (row[j-1] > row[j] and row[j+1] > row[j])):
            print("Found a turning point at coordinates: ({}, {})".format(i, j))
            print("Potential here is {:.3e}".format(V[j,i]))

            ax.scatter(x[j], y[i], V[j,i])


ax.scatter(sep, 0,0, marker='x', color='black')
ax.scatter(0, 0,0, marker='x', color='black')

# ax.contour(x, y, V, Contours, alpha=0.3)
# ax.contourf(x, y, V, L1, colors='red', linestyle='-')

# plt.axis('scaled')
plt.show()
