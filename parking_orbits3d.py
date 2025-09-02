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
import lambertSolver as lt
import planetary_data as pd



# plt.style.use( 'dark_background' )


'''
Setting up Ephemeris Data and Kernels
'''



# Establish Kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")


# Input code - this is variable, please change to whatever you want

planet0 = 'Earth' #Case sensitive
planet1 = 'Jupiter'

#year / month / day
departure0 = '2006-02-01'         # Intial departure date
arrival0 = '2009-01-01'        # Initial arrival date
         
# --------------------



# Frame of reference from the planetary data
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']
mu_earth = 3.986*10**5

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
et_departure = spice.utc2et(departure0)
et_arrival   = spice.utc2et(arrival0)
           
tof = et_arrival - et_departure


'''
Using Lambert's Solver to return the required velocities
'''

          

r_p = 6378.137  # km
h_alt = 1000
i_deg = 80

def circ_velocity_at_point_with_incl(r_point_helio, earth_state_helio, i_deg, mu_central):
    rE = earth_state_helio[:3]
    vE = earth_state_helio[3:]
    
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
        v_helio = earth_state_helio[3:] + v_geo
        

    return v_helio


def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

	
def norm(vec):
    return np.linalg.norm(vec)

earth_centre = calc_ephemeris(departure_planet, et_departure, FRAME, OBSERVER)
    
    
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
x, y, z, theta, phi = sphere_points(earth_centre[:3], radius, num=5)


# Flatten and combine into (N,3) array of points
points = np.vstack((x.flatten(), y.flatten(), z.flatten())).T


# fig = plt.figure(figsize=(8,6))
# ax = fig.add_subplot(111, projection='3d')

# ax.scatter(points[:,0], points[:,1], points[:,2], 
#            c='b', s=20, alpha=0.7)

# # Draw center point
# 

# #ax.set_box_aspect([1,1,1])
# ax.legend()
# plt.show()   
    

states_arrive = calc_ephemeris(arrival_planet, et_arrival, FRAME, OBSERVER)


energies_parking = []
v_parking = np.sqrt(mu_earth / r_p)

length_points = len(points)
points_array = points  # (N,3)
energies_parking = np.zeros(length_points)
v_geo_array = np.zeros((length_points, 3))

for i in range(length_points):
     #v_sc_depart_short = lt.lambert_solver(points[i], states_arrive[:3], tof, sun_mu, trajectory='pro')
     
     
    v_helio_ = circ_velocity_at_point_with_incl(points[i], earth_centre, i_deg, mu_earth)
    
    if v_helio_ is None or np.any(np.isnan(v_helio_)):
        C3_short = np.nan  # or just skip this point


    else:
        v_sc_depart_short = lt.lambert_solver(points[i], states_arrive[:3], tof, sun_mu, trajectory='pro')
        C3_short = min(norm(v_sc_depart_short - v_helio_)**2, cutoff_c3)
        energies_parking[i] = float(C3_short)
        
        v_geo = v_helio_ - earth_centre[3:]
        v_geo_array[i] = v_geo




# vec = points[40] - earth_centre[:3]
# vec_hat = vec / norm(vec)
# print(vec_hat)

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



p = ax.scatter(points_valid[:,0], points_valid[:,1], points_valid[:,2], 
               c=energies_valid, cmap='viridis', s=400)


'''
Arrows to reassure the velocity vector is correct
'''

ax.quiver(
    points_valid[:,0], points_valid[:,1], points_valid[:,2],  # base points
    v_geo_valid[:,0], v_geo_valid[:,1], v_geo_valid[:,2],     # vector components
    length=5000,    # adjust scale for visibility
    color='black',  # or any color
    normalize=True # normalize vectors to the same length
)

# Create sphere mesh
phi, theta = np.linspace(0, 2*np.pi, 30), np.linspace(0, np.pi, 15)
phi, theta = np.meshgrid(phi, theta)
x = earth_centre[0] + r_p * np.sin(theta) * np.cos(phi)
y = earth_centre[1] + r_p * np.sin(theta) * np.sin(phi)
z = earth_centre[2] + r_p * np.cos(theta)

# Plot surface
ax.plot_surface(x, y, z, color='red', alpha=0.1)

# Add colorbar
fig.colorbar(p, ax=ax, shrink=0.6)

ax.set_box_aspect([1,1,1])
ax.set_xlabel('X [km]')
ax.set_ylabel('Y [km]')
ax.set_zlabel('Z [km]')
plt.show()





