# -*- coding: utf-8 -*-
"""
Created on Mon Aug 18 10:30:56 2025

@author: alexa
"""

import numpy as np
import matplotlib.pyplot as plt

#Import Python files

import spiceypy as spice


#Import Lambert tools
import lambertsolve_master as lt
import planetary_data as pd




# Load SPICE kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")

def seconds(days):
    return days * 24 * 3600

def days(seconds_val):
    return seconds_val / (24 * 3600)

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



# Constants
OBSERVER = pd.sun['name']
FRAME = 'ECLIPJ2000'
mu_sun = pd.sun['mu']

# Get planetary data
pd_bodies = pd.bodies
for i in range(10):
    if pd_bodies[i]['name'] == 'Earth':
        pd_earth = pd_bodies[i]
    elif pd_bodies[i]['name'] == 'Jupiter':
        pd_mars = pd_bodies[i]
    elif pd_bodies[i]['name'] == 'Saturn':
        pd_jupiter = pd_bodies[i]

Earth = pd_earth['spice_name']
Jupiter = pd_mars['spice_name']
Saturn = pd_jupiter['spice_name']


#EMJ trajectory times 
time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-10-30"  # Saturn arrival

et_0 = spice.utc2et(time0)
et_1 = spice.utc2et(time1)
et_2 = spice.utc2et(time2)

# LEG 1: Earth to Jupiter

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


print("Gravity assist at Jupiter")
print("------------------------")

print(f"Jupiter V_inf_minis: {vinf_in} km/s")
print(f"Jupiter V_inf_plus: {vinf_out} km/s")
print(f"Deflection angle: {np.degrees(def_angle):.1f} degrees")
# print(f"Heliocentric arriving veloicty: {v_sc_arrive_1} km/s)"
print(f"Heliocentric leaving velocity: {v_sc_depart_2} km/s")
print(f"Leg 1 tof: {days(tof_1):.2f}")
print(f"Leg 2 tof: {days(tof_2):.2f}")

r_jup = 69911
mu_jup = 126.687*10**6

def altitude(delta_a, v_infx):
    e = 1 / np.sin(delta_a / 2)
    r_p = (e - 1) * mu_jup / v_infx**2
    h = r_p - r_jup
    return h

v_inf_mag = norm(vinf_in)
h = altitude(def_angle, v_inf_mag)
print(f"Flyby altitude of Jupiter: {h:.2f} km")


# Saturn arrival analysis
state_Saturn_arrival = calc_ephemeris(Saturn, et_1 + tof_2, FRAME, OBSERVER)
v_inf_Saturn_arrival = v_sc_arrive_2 - state_Saturn_arrival[3:]
    
     





     
























# plt.style.use( 'dark_background' )


'''
Setting up Ephemeris Data and Kernels
'''



# Establish Kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")


# Input code - this is variable, please change to whatever you want

planet0 = 'Jupiter' #Case sensitive
planet1 = 'Saturn'

#year / month / day
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival
         
# --------------------



# Frame of reference from the planetary data
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']


# Loops through the podies to check against the inputted planet0 and planet1
pd_bodies = pd.bodies
for i in range(10):
    if pd_bodies[i]['name'] == planet0:
        pd_req0 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet1:
        pd_req1 = pd_bodies[i]
        
# Holds the spice name of the planets        
departure_planet = pd_req0['spice_name']
arrival_planet = pd_req1['spice_name']



#cut off velocity (if the lambert solver fails to converge)
cutoff_v = 100
cutoff_c3 = cutoff_v ** 2


#Fetch the time data and store into arrays from initial to final
et_departure = spice.utc2et(time1)
et_arrival   = spice.utc2et(time2)
           
tof = et_arrival - et_departure


'''
Using Lambert's Solver to return the required velocities
'''

          

r_p = 69911  # km
h_alt = 576971
i_deg = 80
mu_earth = 3.986*10**5
mu_jup = 126.687*10**6

def circ_velocity_at_point_with_incl(r_point_helio, planet_state_helio, i_deg, mu_central):
    rE = planet_state_helio[:3]
    vE = planet_state_helio[3:]
    
    r_rel = r_point_helio - rE

    r_norm = norm(r_rel)
    r_hat = r_rel / norm(r_rel) 
    
    
    lat_rad = np.arcsin(r_rel[2] / r_norm)  
    
    lat_deg = np.rad2deg(lat_rad)
          
    v_mag = np.sqrt(mu_central / r_norm)
    
    
    i = np.deg2rad(i_deg)
    s_i = np.sin(i)
    c_i = np.cos(i)
    
    val = -r_hat[2] * (c_i / s_i)
    
    
    if abs(lat_deg) > i_deg:
        return np.nan  # infeasible at this point for this i
    else:
        pass

    alpha = np.arctan2(r_hat[1], r_hat[0])
    # two RAAN solutions
    sol = [alpha + np.arcsin(val)]

    
    for Omega in sol:
        # h-hat from (i, Omega)
        h_hat = np.array([np.sin(i)*np.sin(Omega),
                          -np.sin(i)*np.cos(Omega),
                           np.cos(i)])
        t_dir = np.cross(h_hat, r_hat)
        t_norm = np.linalg.norm(t_dir)
        if t_norm < 1e-12:
            continue  # numerical guard
        t_hat = t_dir / t_norm
        
    
        v_geo = v_mag * t_hat         # Earth-centered inertial
        v_helio = planet_state_helio[3:] + v_geo
        

    return v_helio


def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

	
def norm(vec):
    return np.linalg.norm(vec)

planet_centre = calc_ephemeris(departure_planet, et_departure, FRAME, OBSERVER)
    
    
def sphere_points(center, r, num):
    x0, y0, z0 = center[:3]
    theta = np.linspace(0, np.pi, num)      # polar angle
    phi = np.linspace(0, 2*np.pi, num)      # azimuthal angle
    theta, phi = np.meshgrid(theta, phi)

    x = x0 + r * np.sin(theta) * np.cos(phi)
    y = y0 + r * np.sin(theta) * np.sin(phi)
    z = z0 + r * np.cos(theta)

    return x, y, z, theta, phi

radius = r_p + h_alt
x, y, z, theta, phi = sphere_points(planet_centre[:3], radius, num=20)


# Flatten and combine into (N,3) array of points
points = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
states_arrive = calc_ephemeris(arrival_planet, et_arrival, FRAME, OBSERVER)


energies_parking = []
v_parking = np.sqrt(mu_earth / r_p)

length_points = len(points)
points_array = points  # (N,3)
energies_parking = np.zeros(length_points)
v_geo_array = np.zeros((length_points, 3))

for i in range(length_points):
     #v_sc_depart_short = lt.lambert_solver(points[i], states_arrive[:3], tof, sun_mu, trajectory='pro')
     
     
    v_helio_ = circ_velocity_at_point_with_incl(points[i], planet_centre, i_deg, mu_jup)
    
    if v_helio_ is None or np.any(np.isnan(v_helio_)):
        C3_short = np.nan  # or just skip this point


    else:
        v_sc_depart_short = lt.lambert_solver(points[i], states_arrive[:3], tof, sun_mu, trajectory='pro')
        C3_short = min(norm(v_sc_depart_short - v_helio_)**2, cutoff_c3)
        energies_parking[i] = float(C3_short)
        
        v_geo = v_helio_ - planet_centre[3:]
        v_geo_array[i] = v_geo



fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')


# Convert energies_parking to a NumPy array
energies_parking = np.array(energies_parking)
points_array = points  # your (N,3) array of x,y,z

# Create a mask for valid (non-NaN) energies
mask = (~np.isnan(energies_parking)) & (energies_parking > 0) & (~np.isnan(v_geo_array[:,0]))


points_valid = points_array[mask]
energies_valid = energies_parking[mask]
v_geo_valid = v_geo_array[mask]

#32 = 2825.91
#33 = 2867.06

blue_point = np.array([-5.98158e+08,	5.28041e+08,	1.12961e+07])
#33 = -5.98149e+08	5.28048e+08	1.12914e+07



p = ax.scatter(points_valid[:,0], points_valid[:,1], points_valid[:,2], 
               c=energies_valid, cmap='viridis', s=400)


'''
Arrows to reassure the velocity vector is correct
'''

# ax.quiver(
#     points_valid[:,0], points_valid[:,1], points_valid[:,2],  # base points
#     v_geo_valid[:,0], v_geo_valid[:,1], v_geo_valid[:,2],     # vector components
#     length=5000,    # adjust scale for visibility
#     color='black',  # or any color
#     normalize=True # normalize vectors to the same length
# )

vinf_in_hat = vinf_in / np.linalg.norm(vinf_in)


# Create sphere mesh
phi, theta = np.linspace(0, 2*np.pi, 30), np.linspace(0, np.pi, 15)
phi, theta = np.meshgrid(phi, theta)
x = planet_centre[0] + r_p * np.sin(theta) * np.cos(phi)
y = planet_centre[1] + r_p * np.sin(theta) * np.sin(phi)
z = planet_centre[2] + r_p * np.cos(theta)

# Plot surface
ax.plot_surface(x, y, z, color='red', alpha=0.1)

# Add colorbar
fig.colorbar(p, ax=ax, shrink=0.6)

ax.set_box_aspect([1,1,1])

ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
plt.show()








# Inputs
v_inf_in = np.array([-0.33442151, 7.7244149, -1.66190037])  # km/s
v_inf_out = np.array([-7.75599905, -1.40346388, 0.64428571])
r_planet = planet_centre[:3]
mu_jup = 126.687e6
R_jup = 69911
h_min = 10000  # km above surface

v_inf_mag = np.linalg.norm(v_inf_in)
r_p = R_jup + h_min

# Plane normal
n = np.cross(v_inf_in, v_inf_out)
n /= np.linalg.norm(n)

# Periapsis direction
r_hat_p = np.cross(n, v_inf_in)
r_hat_p /= np.linalg.norm(r_hat_p)

# Periapsis position vector in heliocentric frame
r_periapsis = r_planet + r_p * r_hat_p

print("Periapsis position vector (km):", r_periapsis)





