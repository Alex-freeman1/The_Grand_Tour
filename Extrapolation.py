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
planet1 = 'Jupiter'
planet2 = 'Saturn'
planet3 = 'URANUS BARYCENTER'
planet4 = 'NEPTUNE BARYCENTER'

orbital_periods = {
    'EARTH': 365.25 * 86400,  # seconds in one Earth year
    'JUPITER': 4333 * 86400,
    'SATURN': 10759 * 86400,
    'URANUS': 30687 * 86400,
    'NEPTUNE': 60190 * 86400
}

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
FRAME='ECLIPJ2000'
mu_sun = pd.sun['mu']

def norm(vec):
    return np.linalg.norm(vec)

def radians(deg):
    rads = deg * (pi/180)
    return rads

def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

# Year - Month - Day
departure0 = '1977-08-20'         
arrival0 = '1979-07-09' 
arrival1 = '1980-11-12'
arrival2 = '1986-01-24' 
arrival3 = '1990-08-25'         

et_0 = spice.utc2et(departure0)
et_1 = spice.utc2et(arrival0)
et_2 = spice.utc2et(arrival1)
et_3 = spice.utc2et(arrival2)
et_4 =  spice.utc2et(arrival3)

et_times = [et_0, et_1, et_2, et_3, et_4]

#string_var = "-705844751.8	-657028750.8 	-602686124.7"

string_var = "-705844751.8 -661348750.8 -611537698.9"

arr = np.array(string_var.split(), dtype=float)
et_0 = arr[0]
et_1 = arr[1]
et_2 = arr[2]



ephem_0 = calc_ephemeris(departure_planet, et_0, FRAME, OBSERVER)
ephem_1 = calc_ephemeris(arrival_planet0, et_1, FRAME, OBSERVER)
ephem_2 = calc_ephemeris(arrival_planet1, et_2, FRAME, OBSERVER)
ephem_3 = calc_ephemeris('URANUS BARYCENTER', et_times[3], FRAME, OBSERVER)
ephem_4 = calc_ephemeris('NEPTUNE BARYCENTER', et_times[4], FRAME, OBSERVER)

tofs = []

for i in range(len(et_times) - 1):
    tofs.append(et_times[i+1] - et_times[i])


def lambert_sol(dep, arr, tof):
    try:
        v_sc_depart_short, v_sc_arrive_short = lt.lambert_solver(
            # states_arrive is a 6 element vector, the first three represent the position vector: x,y,z
            # the final three represent the velocity vector: x,y,z.
            # therefore :3 is taking the first three position vectors passing it as R1
            dep[ :3 ],
            arr[ :3 ], #R2
            tof,
            mu_sun,
            trajectory='pro'
        )
    except Exception as e:
        print(f"Lambert solution failed: {e}")
        raise
    
    return v_sc_depart_short



def two_body(t, y):
    r = y[:3]
    v = y[3:]
    norm_r = np.linalg.norm(r)
    a = -mu_sun * r / norm_r**3
    return np.concatenate((v, a))


ephems = [ephem_0, ephem_1, ephem_2, ephem_3, ephem_4]
n_legs = len(ephems) - 1  # 3 legs
X_list = []  # Initial state vectors
t_spans = []
t_evals = []
solutions = []

# Generate initial state vectors and time spans
for i in range(n_legs):
    r = ephems[i][:3]
    v = lambert_sol(ephems[i], ephems[i+1], tofs[i])
    X = list(r) + list(v)
    X_list.append(X)
    
    dt = et_times[i+1] - et_times[i]
    t_spans.append((0, dt))
    t_evals.append(np.linspace(0, dt, 1000))

# Solve using solve_ivp
for i in range(n_legs):
    sol = solve_ivp(
        two_body, t_spans[i], X_list[i], t_eval=t_evals[i],
        method='RK45', rtol=1e-12, atol=1e-12)
    solutions.append(sol)
    
solution0, solution1, solution2, solution3 = solutions

step = 10000  # in seconds
et0 = spice.utc2et(departure0)  # reference start time


et_earth = np.arange(et0, et0 + orbital_periods['EARTH'], step)
et_jupiter  = np.arange(et0, et0 + orbital_periods['JUPITER'],  step)
et_saturn = np.arange(et0, et0 + orbital_periods['SATURN'],  step)
et_uranus = np.arange(et0, et0 + orbital_periods['URANUS'],  step)
et_neptune = np.arange(et0, et0 + orbital_periods['NEPTUNE'],  step)



r_earth = np.array([calc_ephemeris(departure_planet, t, FRAME, OBSERVER)[:3] for t in et_earth])
r_jupiter  = np.array([calc_ephemeris('JUPITER BARYCENTER',  t, FRAME, OBSERVER)[:3] for t in et_jupiter])
r_saturn = np.array([calc_ephemeris('SATURN BARYCENTER',  t, FRAME, OBSERVER)[:3] for t in et_saturn])
r_uranus = np.array([calc_ephemeris('URANUS BARYCENTER',  t, FRAME, OBSERVER)[:3] for t in et_uranus])
r_neptune = np.array([calc_ephemeris('NEPTUNE BARYCENTER',  t, FRAME, OBSERVER)[:3] for t in et_neptune])


# Extract x, y, z components
x_earth, y_earth, z_earth = r_earth[:, 0], r_earth[:, 1], r_earth[:, 2]
x_jupiter, y_jupiter, z_jupiter = r_jupiter[:, 0], r_jupiter[:, 1], r_jupiter[:, 2]
x_saturn, y_saturn, z_saturn = r_saturn[:, 0], r_saturn[:, 1], r_saturn[:, 2]
x_uranus, y_uranus, z_uranus = r_uranus[:, 0], r_uranus[:, 1], r_uranus[:, 2]
x_neptune, y_neptune, z_neptune = r_neptune[:, 0], r_neptune[:, 1], r_neptune[:, 2]

# Extract components
x, y, z = solution0.y[0], solution0.y[1], solution0.y[2]
xj1, yj1, zj1 = solution1.y[0], solution1.y[1], solution1.y[2]
xs1, ys2, zs3 = solution2.y[0], solution2.y[1], solution2.y[2]
xu1, yu2, zu3 = solution3.y[0], solution3.y[1], solution3.y[2]




fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

max_range = np.max(np.abs(x_jupiter))
ax.set_zlim([-max_range, max_range])


r0_earth = ephems[0][:3] 
r1_jupiter = ephems[1][:3]
r2_saturn = ephems[2][:3]
r3_uranus = ephems[3][:3]
r4_neptune = ephems[4][:3]

# Plot planet's orbit
ax.plot(x_earth, y_earth, z_earth, label="Earth's Orbit", color="b")
ax.plot(x_jupiter, y_jupiter, z_jupiter, label="Jupiter's Orbit", color="r")
ax.plot(x_saturn, y_saturn, z_saturn, label="Saturn's Orbit", color="g")
# ax.plot(x_uranus, y_uranus, z_uranus, label="Uranus's Orbit", color="c")
# ax.plot(x_neptune, y_neptune, z_neptune, label="Neptunes's Orbit", color="b")

#Mark the Sun at (0,0,0)# Mark the Sun at (0,0,0)
ax.scatter(0, 0, 0, color='yellow', s=100, label="Sun")


ax.plot(x, y,z, label='Lambert Trajectory')
ax.plot(xj1, yj1, zj1, label='Lambert Trajectory')
# ax.plot(xs1, ys2, zs3, label='Lambert Trajectory')
# ax.plot(xu1, yu2, zu3, label='Lambert Trajectory')

ax.scatter(*r0_earth, color='blue', marker='o', s=100, label='Earth Departure')
ax.scatter(*r1_jupiter, color='red', marker='o', s=100, label='Jupiter Arrival')
ax.scatter(*r2_saturn, color='green', marker='o', s=100, label='Saturn Arrival')
# ax.scatter(*r3_uranus, color='pink', marker='o', s=100, label='Uranus Arrival')
# ax.scatter(*r4_neptune, color='blue', marker='o', s=100, label='Neptune Arrival')
ax.set_xlabel("X (km)")
ax.set_ylabel("Y (km)")
ax.set_zlabel("Z (km)")

ax.legend(fontsize= 25)

plt.show()


# ------------------------


