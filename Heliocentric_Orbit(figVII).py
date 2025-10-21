
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 17:07:55 2025
@author: alexa
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import spiceypy as spice
import lambertsolve_master as lt
import planetary_data as pd
spice.furnsh("data\de440.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")
planet0 = 'Earth' #Case sensitive
planet1 = 'Mars'
planet2 = 'Jupiter'
pd_bodies = pd.bodies
for i in range(10):
    if pd_bodies[i]['name'] == planet0:
        pd_req0 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet1:
        pd_req1 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet2:
        pd_req2 = pd_bodies[i]
        
# Holds the spice name of the planets        
departure_planet = pd_req0['spice_name']
arrival_planet0 = pd_req1['spice_name']
arrival_planet1 = pd_req2['spice_name']
pi = np.pi
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000 '
mu_sun = pd.sun['mu']
def norm(vec):
    return np.linalg.norm(vec)
def radians(deg):
    rads = deg * (pi/180)
    return rads
def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])
# Year - Month - Day
departure0 = '1965-01-26'          # Bad transfer
arrival0 = '1965-12-08'  
arrival1 = '1966-08-23'    

# departure0 = '1965-01-31'         
# arrival0 = '1966-03-08'  
# arrival1 = '1969-02-24'  
   
et_departure = spice.utc2et(departure0)
et_arrival = spice.utc2et(arrival0)
et_arrival1 = spice.utc2et(arrival1)
ephem_departure = calc_ephemeris(departure_planet, et_departure, FRAME, OBSERVER)
ephem_arrival = calc_ephemeris(arrival_planet0, et_arrival, FRAME, OBSERVER)
ephem_arrival1 = calc_ephemeris(arrival_planet1, et_arrival1, FRAME, OBSERVER)
tof1 = et_arrival - et_departure
tof2 = et_arrival1 - et_arrival
try:
    v_sc_depart_short, v_sc_arrive_short = lt.lambert_solver(
        # states_arrive is a 6 element vector, the first three represent the position vector: x,y,z
        # the final three represent the velocity vector: x,y,z.
        # therefore :3 is taking the first three position vectors passing it as R1
        ephem_departure[ :3 ],
        ephem_arrival[ :3 ], #R2
        tof1,
        mu_sun,
        trajectory='pro'
    )
except:
    v_sc_depart_short = np.array( [1000, 1000, 1000] )
    v_sc_arrive_short = np.array( [1000, 1000, 1000] )
 
x0, y0, z0 = ephem_departure[:3] 
vx0, vy0, vz0 = v_sc_depart_short
x1, y1, z1 = ephem_arrival[:3] 
try:
    v_sc_depart_short, v_sc_arrive_short = lt.lambert_solver(
        # states_arrive is a 6 element vector, the first three represent the position vector: x,y,z
        # the final three represent the velocity vector: x,y,z.
        # therefore :3 is taking the first three position vectors passing it as R1
        ephem_arrival[ :3 ],
        ephem_arrival1[ :3 ], #R2
        tof2,
        mu_sun,
        trajectory='pro'
    )
except:
    v_sc_depart_short = np.array( [1000, 1000, 1000] )
    v_sc_arrive_short = np.array( [1000, 1000, 1000] )
vx1, vy1, vz1 = v_sc_depart_short
x2,y2,z2 = ephem_arrival1[:3] 
# Initial state vector
X0 = [x0, y0, z0, vx0, vy0, vz0]   
X1 = [x1, y1, z1, vx1, vy1, vz1]
def two_body(t, y):
    r = y[:3]
    v = y[3:]
    norm_r = np.linalg.norm(r)
    a = -mu_sun * r / norm_r**3
    return np.concatenate((v, a))
t_span0 = (0, (et_arrival - et_departure))  # assuming t1, t2 are in seconds
t_eval0 = np.linspace(t_span0[0], t_span0[1], 1000) 
t_span1 = (0, (et_arrival1 - et_arrival))  # assuming t1, t2 are in seconds
t_eval1 = np.linspace(t_span1[0], t_span1[1], 1000) 
# Solve the system using solve_ivp with specified tolerances using RK45
solution0 = solve_ivp(
    two_body, t_span0, X0, t_eval=t_eval0, method='RK45',
    rtol=1e-12, atol=1e-12)
solution1 = solve_ivp(
    two_body, t_span1, X1, t_eval=t_eval1, method='RK45',
    rtol=1e-12, atol=1e-12)
step = 1000  # in seconds
et0 = spice.utc2et(departure0)  # reference start time
orbital_periods = {
    'EARTH': 365.25 * 86400,  # seconds in one Earth year
    'MARS': 686.98 * 86400,   # seconds in one Mars year
    'JUPITER': 4333 * 86400
}
et_earth = np.arange(et0, et0 + orbital_periods['EARTH'], step)
et_mars  = np.arange(et0, et0 + orbital_periods['MARS'],  step)
et_jupiter = np.arange(et0, et0 + orbital_periods['JUPITER'],  step)
r_earth = np.array([calc_ephemeris(departure_planet, t, FRAME, OBSERVER)[:3] for t in et_earth])
r_mars  = np.array([calc_ephemeris(arrival_planet0,  t, FRAME, OBSERVER)[:3] for t in et_mars])
r_jupiter = np.array([calc_ephemeris('JUPITER BARYCENTER',  t, FRAME, OBSERVER)[:3] for t in et_jupiter])
# Extract x, y, z components
x_earth, y_earth, z_earth = r_earth[:, 0], r_earth[:, 1], r_earth[:, 2]
x_mars, y_mars, z_mars = r_mars[:, 0], r_mars[:, 1], r_mars[:, 2]
x_jupiter, y_jupiter, z_jupiter = r_jupiter[:, 0], r_jupiter[:, 1], r_jupiter[:, 2]
# Extract components
x, y, z = solution0.y[0], solution0.y[1], solution0.y[2]
xj1, yj1, zj1 = solution1.y[0], solution1.y[1], solution1.y[2]
#vx, vy, vz = solution.y[3], solution.y[4], solution.y[5]
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')
max_range = np.max(np.abs(x_jupiter))
ax.set_zlim([-max_range, max_range])
r0_earth = x0, y0, z0 
r1_mars = x1,y1,z1
r2_jupiter = x2,y2,z2 
# Plot planet's orbit
ax.plot(x_earth, y_earth, z_earth, label="Earth's Orbit", color="b")
ax.plot(x_mars, y_mars, z_mars, label="Mar's Orbit", color="r")
ax.plot(x_jupiter, y_jupiter, z_jupiter, label="Jupiter's Orbit", color="g")
#Mark the Sun at (0,0,0)# Mark the Sun at (0,0,0)
ax.scatter(0, 0, 0, color='yellow', s=100, label="Sun")
ax.plot(x, y,z, label='Lambert Trajectory')
ax.plot(xj1, yj1, zj1, label='Lambert Trajectory')
ax.scatter(*r0_earth, color='blue', marker='o', s=100, label='Earth Departure')
ax.scatter(*r1_mars, color='red', marker='o', s=100, label='Mars Arrival')
ax.scatter(*r2_jupiter, color='green', marker='o', s=100, label='Jupiter Arrival')
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_zlabel("Z (km)")
ax.legend(fontsize = 25)
plt.show()
# ------------------------