# -*- coding: utf-8 -*-
"""
Created on Thu Jun 26 13:43:40 2025

@author: alexa
"""

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import lambertsolve_master as lt
import planetary_data as pd





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

def norm(vec):
    return np.linalg.norm(vec)

def calc_vinfinity(tof, target_planet, departure_planet, et0, vinf_target):
    """
    Function to compute difference between target vinfinity magnitude and
    current vinfinity magnitude for a given time of flight (tof).
    """
    r1_target = calc_ephemeris(target_planet, et0 + tof, FRAME, OBSERVER)[:3]
    state_departure = calc_ephemeris(departure_planet, et0, FRAME, OBSERVER)
    
    v0_sc_depart, v1_sc_arrive = lt.lambert_solver(
        state_departure[:3],    
        r1_target,   
        tof,
        mu_sun,
        trajectory='pro'
    )
    
    vinf = np.linalg.norm(v0_sc_depart - state_departure[3:])
    return vinf_target - vinf

def vinfinity_match(planet0, planet1, v_sc_incoming, et0, tof0, diff_step=1e-3, tol=1e-4):
    """
    Match V-infinity magnitudes for hyperbolic flyby
    """
    state0_planet0 = calc_ephemeris(planet0, et0, FRAME, OBSERVER)
    vinf_incoming = v_sc_incoming - state0_planet0[3:]
    vinf_magnitude = np.linalg.norm(vinf_incoming)
    
   
    def root_func(tof):
        return calc_vinfinity(tof, planet1, planet0, et0, vinf_magnitude)
    
    tof, steps = lt.newton_root_single_fd(root_func, tof0, tol=tol, diff_step=diff_step)
    
    r1_planet1 = calc_ephemeris(planet1, et0 + tof, FRAME, OBSERVER)[:3]
    
    v_sc_depart, v_sc_arrive = lt.lambert_solver(
        state0_planet0[:3],    
        r1_planet1,   
        tof,
        mu_sun,
        trajectory='pro'
    )
    
    return tof, v_sc_depart, v_sc_arrive

# time0 = "1977-08-19"  # Earth departure
# time1 = "1979-01-16"  # Jupiter flyby
# time2 = "1980-08-15"  # Saturn arrival

# et_0 = spice.utc2et(time0)
# et_1 = spice.utc2et(time1)
# et_2 = spice.utc2et(time2)

#print(et_0, et_1, et_2)

string_var = "-705844751.8 -661348750.8 -611537698.9"

arr = np.array(string_var.split(), dtype=float)
et_0 = arr[0]
et_1 = arr[1]
et_2 = arr[2]



state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
state1 = calc_ephemeris(Jupiter, et_1, FRAME, OBSERVER)
tof_1 = et_1 - et_0



# Lambert solution for Earth-Jupiter leg
v_sc_depart_1, v_sc_arrive_1 = lt.lambert_solver(
    state0[:3],    
    state1[:3],   
    tof_1,
    mu_sun,
    trajectory='pro'
)

#  ----------------------------------------------------------------------------

# LEG 2: Jupiter to Saturn 
tof_guess_2 = et_2 - et_1  # Initial guess


tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
    Jupiter, Saturn,
    v_sc_arrive_1,  # spacecraft velocity at Jupiter arrival
    et_1, tof_guess_2
)


# Verify the Jupiter flyby
state_Jupiter_flyby = calc_ephemeris(Jupiter, et_1, FRAME, OBSERVER)
vinf_in = v_sc_arrive_1 - state_Jupiter_flyby[3:]
vinf_out = v_sc_depart_2 - state_Jupiter_flyby[3:]

# Calculate deflection angle
def_angle = np.arccos(np.dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))



mu_jup = 1.2668653e8  # km^3/s^2
r_jup = 71492.0  # km



def altitude(delta_a, v_infx):
    e = 1 / np.sin(delta_a / 2)
    r_p = (e - 1) * mu_jup / v_infx**2
    h = r_p - r_jup
    return h


v_inf_mag = norm(vinf_in)
h = altitude(def_angle, v_inf_mag)
print(h)
print(h + r_jup)
if h < 0:
    raise ValueError('h should not be negative')

rp = r_jup + h
v_inf_mag = np.linalg.norm(vinf_in)
a = -mu_jup / v_inf_mag**2


e = 1 + rp * v_inf_mag**2 / mu_jup
delta = 2 * np.arcsin(1 / e)
print(np.rad2deg(delta))
def rotate_vector(v, axis, angle):
    axis = axis / np.linalg.norm(axis)
    return (v * np.cos(angle) +
            np.cross(axis, v) * np.sin(angle) +
            axis * np.dot(axis, v) * (1 - np.cos(angle)))

# Arbitrary impact parameter vector perpendicular to v_inf_in
impact_axis = np.cross(vinf_in, [0, 0, 1])
if np.linalg.norm(impact_axis) < 1e-6:
    impact_axis = np.cross(vinf_in, [0, 1, 0])

impact_axis = impact_axis / np.linalg.norm(impact_axis)

v_inf_out = rotate_vector(vinf_in, impact_axis, delta)


# Generate true anomalies
theta = np.linspace(-np.radians(120), np.radians(120), 1000)  

# Hyperbolic trajectory in perifocal frame (orbital plane)
r = a * (e**2 - 1) / (1 + e * np.cos(theta))
x_pf = r * np.cos(theta)
y_pf = r * np.sin(theta)
z_pf = np.zeros_like(x_pf)

# Stack into array
trajectory_pf = np.vstack([x_pf, y_pf, z_pf])  # 3 x N

# Define perifocal axes
x_hat = vinf_in / np.linalg.norm(vinf_in)
z_vec = np.cross(vinf_in, v_inf_out)
z_hat = z_vec / np.linalg.norm(z_vec)
y_vec = np.cross(z_hat, x_hat)
y_hat = y_vec / np.linalg.norm(y_vec)

# Rotation matrix: columns are basis vectors of perifocal frame
R = np.column_stack((x_hat, y_hat, z_hat))

# Transform to inertial frame
trajectory_inertial = R @ trajectory_pf

def align_plane(trajectory, plane_normal):
   
    target_normal = np.array([0, 0, 1.0])

    # Compute rotation axis (cross product)
    axis = np.cross(plane_normal, target_normal)
    if np.linalg.norm(axis) < 1e-6:
        return trajectory  

    axis /= np.linalg.norm(axis)
    angle = np.arccos(np.dot(plane_normal, target_normal))

    rotated_traj = np.array([rotate_vector(p, axis, angle) for p in trajectory.T]).T
    return rotated_traj

# Apply this right before plotting:
trajectory_aligned = align_plane(trajectory_inertial, z_hat)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot(*trajectory_aligned, color='blue', label='Hyperbolic Flyby')
# ax.quiver(0, 0, 0, *vinf_in * r_jup, color='green', label='v_inf_in')
# ax.quiver(0, 0, 0, *v_inf_out * r_jup, color='red', label='v_inf_out')


# Draw Jupiter as a sphere
u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
xj = r_jup * np.cos(u) * np.sin(v)
yj = r_jup * np.sin(u) * np.sin(v)
zj = r_jup * np.cos(v)
ax.plot_surface(xj, yj, zj, color='orange', alpha=0.5)

# Set equal aspect ratio (spherical appearance)
max_range = r_jup * 10# enough to include the hyperbolic arc
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