# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 13:43:40 2025

@author: alexa
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import lambertsolve_master as lt
import planetary_data as pd
from scipy.integrate import solve_ivp




# Load SPICE kernels
spice.furnsh("data/de440.bsp")
spice.furnsh("data/latest_leapseconds (1).tls")

# Planet setup
planet0 = 'Earth'
planet1 = 'Jupiter'
planet2 = 'Saturn'
pd_bodies = pd.bodies

for i in range(10):
    if pd_bodies[i]['name'] == planet0:
        pd_req0 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet1:
        pd_req1 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet2:
        pd_req2 = pd_bodies[i]
        
# Holds the spice name of the planets        
Earth = pd_req0['spice_name']
Jupiter = pd_req1['spice_name']
Saturn = pd_req2['spice_name']



pi = np.pi
OBSERVER = pd.sun['name']
FRAME = 'ECLIPJ2000'
mu_sun = pd.sun['mu']


def norm(vec):
    return np.linalg.norm(vec)

def radians(deg):
    rads = deg * (pi/180)
    return rads

AU = 1.496e+8
def seconds(days):
    return days * 24 * 3600

def calc_ephemeris(target, ets, frame, observer):
    return np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])

time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival


vinf_in = [-0.33442151,  7.7244149 , -1.66190037]
vinf_out= [-7.75599905, -1.40346388,  0.64428571]
v_sc_arrive_1 = [-9.14381453, -1.47503865, -1.426696]
v_sc_depart_2 = [-16.56539206, -10.60291743,  0.87949009]

v_inf_in = np.array([-0.33442151,  7.7244149 , -1.66190037]) 
h = 1901522  # altitude in km above Jupiter
mu_J = 1.2668653e8  # km^3/s^2
R_J = 71492.0  # km
rp = R_J + h


v_inf_mag = np.linalg.norm(v_inf_in)
a = -mu_J / v_inf_mag**2


e = 1 + rp * v_inf_mag**2 / mu_J
delta = 2 * np.arcsin(1 / e)

def rotate_vector(v, axis, angle):
    axis = axis / np.linalg.norm(axis)
    return (v * np.cos(angle) +
            np.cross(axis, v) * np.sin(angle) +
            axis * np.dot(axis, v) * (1 - np.cos(angle)))

# Arbitrary impact parameter vector perpendicular to v_inf_in
impact_axis = np.cross(v_inf_in, [0, 0, 1])
if np.linalg.norm(impact_axis) < 1e-6:
    impact_axis = np.cross(v_inf_in, [0, 1, 0])  # handle edge case

impact_axis = impact_axis / np.linalg.norm(impact_axis)

v_inf_out = rotate_vector(v_inf_in, impact_axis, delta)


# Generate true anomalies
theta = np.linspace(-np.radians(80), np.radians(80), 500)  # near-perijove arc

# Hyperbolic trajectory in perifocal frame (orbital plane)
r = a * (e**2 - 1) / (1 + e * np.cos(theta))
x_pf = r * np.cos(theta)
y_pf = r * np.sin(theta)
z_pf = np.zeros_like(x_pf)

# Stack into array
trajectory_pf = np.vstack([x_pf, y_pf, z_pf])  # 3 x N

# Define perifocal axes
x_hat = v_inf_in / np.linalg.norm(v_inf_in)
z_hat = np.cross(v_inf_in, v_inf_out)
z_hat = z_hat / np.linalg.norm(z_hat)
y_hat = np.cross(z_hat, x_hat)

# Rotation matrix: columns are basis vectors of perifocal frame
R = np.column_stack((x_hat, y_hat, z_hat))

# Transform to inertial frame
trajectory_inertial = R @ trajectory_pf


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(*trajectory_inertial, color='blue', label='Hyperbolic Flyby')
# ax.quiver(0, 0, 0, *v_inf_in * R_J, color='green', label='v_inf_in')
# ax.quiver(0, 0, 0, *v_inf_out * R_J, color='red', label='v_inf_out')

# Draw Jupiter as a sphere
u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
xj = R_J * np.cos(u) * np.sin(v)
yj = R_J * np.sin(u) * np.sin(v)
zj = R_J * np.cos(v)
ax.plot_surface(xj, yj, zj, color='orange', alpha=0.5)

# Set equal aspect ratio (spherical appearance)
max_range = R_J * 30  # enough to include the hyperbolic arc
ax.set_xlim(-max_range, max_range)
ax.set_ylim(-max_range, max_range)
ax.set_zlim(-max_range, max_range)
ax.set_box_aspect([1, 1, 1])  # Equal scale for all axes

ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
ax.legend()
ax.set_box_aspect([1,1,1])
plt.title('Gravity Assist at Jupiter (to Saturn)')
plt.show()