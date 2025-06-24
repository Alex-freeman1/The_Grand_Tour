# -*- coding: utf-8 -*-
"""
Created on Thu May 29 23:29:40 2025

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
Earth = pd_req0['spice_name']
Mars = pd_req1['spice_name']
Jupiter = pd_req2['spice_name']

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

# Mission parameters
departure0 = '2006-02-15'  # Earth departure
arrival0 = '2007-01-10'    # Mars arrival
arrival1 = '2010-01-01'    # Jupiter arrival

et_departure = spice.utc2et(departure0)
et_arrival = spice.utc2et(arrival0)
et_arrival1 = spice.utc2et(arrival1)

tof1 = et_arrival - et_departure      # Earth to Mars time
tof2 = et_arrival1 - et_arrival       # Mars to Jupiter time

print(f"Earth to Mars: {tof1/86400:.1f} days")
print(f"Mars to Jupiter: {tof2/86400:.1f} days")

# Get ephemeris data for all mission points
ephem_departure = calc_ephemeris(Earth, et_departure, FRAME, OBSERVER)
ephem_arrival = calc_ephemeris(Mars, et_arrival, FRAME, OBSERVER)
ephem_arrival1 = calc_ephemeris(Jupiter, et_arrival1, FRAME, OBSERVER)

# Solve Lambert's problem for Earth to Mars
try:
    v_sc_depart, v_sc_arrive = lt.lambert_solver(
        ephem_departure[:3],  # Earth position at departure
        ephem_arrival[:3],    # Mars position at arrival
        tof1,                 # Time of flight in seconds
        mu_sun,               # Sun's gravitational parameter
        trajectory='retro'      # Prograde trajectory
    )
    print("Earth-Mars Lambert solver successful")
except Exception as e:
    print(f"Earth-Mars Lambert solver failed: {e}")
    v_sc_depart = np.array([1000, 1000, 1000])
    v_sc_arrive = np.array([1000, 1000, 1000])

# Solve Lambert's problem for Mars to Jupiter
try:
    v_sc_depart2, v_sc_arrive2 = lt.lambert_solver(
        ephem_arrival[:3],    # Mars position at departure
        ephem_arrival1[:3],   # Jupiter position at arrival
        tof2,                 # Time of flight in seconds
        mu_sun,               # Sun's gravitational parameter
        trajectory='retro'      # Prograde trajectory
    )
    print("Mars-Jupiter Lambert solver successful")
except Exception as e:
    print(f"Mars-Jupiter Lambert solver failed: {e}")
    v_sc_depart2 = np.array([1000, 1000, 1000])
    v_sc_arrive2 = np.array([1000, 1000, 1000])

# Initial conditions for both trajectory segments
X0 = [*ephem_departure[:3], *v_sc_depart]      # Earth to Mars
X1 = [*ephem_arrival[:3], *v_sc_depart2]       # Mars to Jupiter

def two_body(t, y):
    """Two-body dynamics for spacecraft motion"""
    r = y[:3]
    v = y[3:]
    norm_r = np.linalg.norm(r)
    a = -mu_sun * r / norm_r**3
    return np.concatenate((v, a))

# Integrate first trajectory (Earth to Mars)
t_span1 = (0, tof1)
t_eval1 = np.linspace(0, tof1, 1000)
solution1 = solve_ivp(two_body, t_span1, X0, t_eval=t_eval1, method='RK45', rtol=1e-12, atol=1e-12)

# Integrate second trajectory (Mars to Jupiter)
t_span2 = (0, tof2)
t_eval2 = np.linspace(0, tof2, 1000)
solution2 = solve_ivp(two_body, t_span2, X1, t_eval=t_eval2, method='RK45', rtol=1e-12, atol=1e-12)

# Extract trajectory data
if solution1.success:
    x_traj1, y_traj1, z_traj1 = solution1.y[0], solution1.y[1], solution1.y[2]
    print("First trajectory integration successful")
else:
    x_traj1 = y_traj1 = z_traj1 = np.array([0])
    print("First trajectory integration failed")

if solution2.success:
    x_traj2, y_traj2, z_traj2 = solution2.y[0], solution2.y[1], solution2.y[2]
    print("Second trajectory integration successful")
else:
    x_traj2 = y_traj2 = z_traj2 = np.array([0])
    print("Second trajectory integration failed")


spacecraft_pos1 = np.stack((x_traj1, y_traj1, z_traj1), axis=-1)
spacecraft_pos2 = np.stack((x_traj2, y_traj2, z_traj2), axis=-1)


# Setup orbital data for planets
start = '2005-11-26' 
et_start = int(spice.utc2et(start))

# Use a common simulation period (longer to show full mission)
#simulation_period = int(seconds(365.25) * 6)  # 6 years
time_step = 86400  # 1 day

# Calculate transfer timing
transfer1_start_day = int((et_departure - et_start) / time_step)
transfer1_duration_days = int(tof1 / time_step)
transfer2_start_day = int((et_arrival - et_start) / time_step)
transfer2_duration_days = int(tof2 / time_step)

# Generate planet positions for the same time period
mars_pos = []
earth_pos = []
jupiter_pos = []


earth_orbital_period = seconds(365.25)
mars_orbital_period = seconds(687)
juptier_orbital_period = seconds(4333)

et_end = et_start + int(earth_orbital_period) 
et_end_mars = et_start + int(mars_orbital_period) 
et_end_jupiter = et_start + int(juptier_orbital_period) 

for t in range(et_start, et_end, time_step):  # Use time_step instead of every second
    earth_pos.append(calc_ephemeris(Earth, t, FRAME, OBSERVER)[:3])
    
for t1 in range(et_start, et_end_mars, time_step):
    mars_pos.append(calc_ephemeris(Mars, t1, FRAME, OBSERVER)[:3])

for t2 in range(et_start, et_end_jupiter, time_step):
    jupiter_pos.append(calc_ephemeris(Jupiter, t2, FRAME, OBSERVER)[:3])

mars_pos = np.array(mars_pos)
earth_pos = np.array(earth_pos)
jupiter_pos = np.array(jupiter_pos)



# Set up the 3D figure and axis
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# Set axis limits to include Jupiter
margin = 1.2
max_range = max(np.max(np.abs(mars_pos)), np.max(np.abs(earth_pos)), np.max(np.abs(jupiter_pos)))
ax.set_xlim(-max_range * margin, max_range * margin)
ax.set_ylim(-max_range * margin, max_range * margin)
ax.set_zlim(-max_range * margin, max_range * margin)

# Labels and title
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('Earth-Mars-Jupiter Transfer Trajectory')

# Create plot objects
earth_dot, = ax.plot([], [], [], 'bo', markersize=8, label='Earth')
mars_dot, = ax.plot([], [], [], 'ro', markersize=8, label='Mars')
jupiter_dot, = ax.plot([], [], [], 'orange', marker='o', markersize=10, label='Jupiter')
spacecraft_dot, = ax.plot([], [], [], 'go', markersize=6, label='Spacecraft')
sun_dot, = ax.plot([0], [0], [0], 'yo', markersize=12, label='Sun')

# Trail objects - FIXED: Added comma for jupiter_trail
earth_trail, = ax.plot([], [], [], 'b-', alpha=0.3, linewidth=1)
mars_trail, = ax.plot([], [], [], 'r-', alpha=0.3, linewidth=1)
jupiter_trail, = ax.plot([], [], [], 'orange', alpha=0.3, linewidth=1)  # FIXED
transfer_trail1, = ax.plot([], [], [], 'g-', alpha=0.8, linewidth=2, label='Earth-Mars')
transfer_trail2, = ax.plot([], [], [], 'm-', alpha=0.8, linewidth=2, label='Mars-Jupiter')

# Trail data
trail_length = 50
earth_trail_data = [[], [], []]
mars_trail_data = [[], [], []]
jupiter_trail_data = [[], [], []]
transfer_trail1_data = [[], [], []]
transfer_trail2_data = [[], [], []]

def animate(frame):
    """Animation function"""
    # Handle wraparound for long simulations
    earth_idx = frame % len(earth_pos)
    mars_idx = frame % len(mars_pos)
    jupiter_idx = frame % len(jupiter_pos)
    
    # Current planet positions
    earth_current = earth_pos[earth_idx]
    mars_current = mars_pos[mars_idx]
    jupiter_current = jupiter_pos[jupiter_idx]
    
    # Update planet positions - FIXED: Jupiter position update
    earth_dot.set_data_3d([earth_current[0]], [earth_current[1]], [earth_current[2]])
    mars_dot.set_data_3d([mars_current[0]], [mars_current[1]], [mars_current[2]])
    jupiter_dot.set_data_3d([jupiter_current[0]], [jupiter_current[1]], [jupiter_current[2]])  # FIXED
    
    if frame % 3 == 0:
       for i, pos in enumerate(earth_current):
           earth_trail_data[i].append(pos)
           if len(earth_trail_data[i]) > trail_length:
               earth_trail_data[i].pop(0)

       for i, pos in enumerate(mars_current):
           mars_trail_data[i].append(pos)
           if len(mars_trail_data[i]) > trail_length:
               mars_trail_data[i].pop(0)

       for i, pos in enumerate(jupiter_current):
           jupiter_trail_data[i].append(pos)
           if len(jupiter_trail_data[i]) > trail_length:
               jupiter_trail_data[i].pop(0)

       earth_trail.set_data_3d(earth_trail_data[0], earth_trail_data[1], earth_trail_data[2])
       mars_trail.set_data_3d(mars_trail_data[0], mars_trail_data[1], mars_trail_data[2])
       jupiter_trail.set_data_3d(jupiter_trail_data[0], jupiter_trail_data[1], jupiter_trail_data[2])

    
    earth_trail.set_data_3d(earth_trail_data[0], earth_trail_data[1], earth_trail_data[2])
    mars_trail.set_data_3d(mars_trail_data[0], mars_trail_data[1], mars_trail_data[2])
    jupiter_trail.set_data_3d(jupiter_trail_data[0], jupiter_trail_data[1], jupiter_trail_data[2])  # FIXED: Added this line
    
    # Handle spacecraft trajectory during transfer windows
    spacecraft_visible = False
    
    # First transfer: Earth to Mars
    if transfer1_start_day <= frame < transfer1_start_day + transfer1_duration_days:
        transfer_frame = frame - transfer1_start_day
        traj_idx = int((transfer_frame / transfer1_duration_days) * len(x_traj1))
        traj_idx = min(traj_idx, len(x_traj1) - 1)
        
        spacecraft_pos = [x_traj1[traj_idx], y_traj1[traj_idx], z_traj1[traj_idx]]
        spacecraft_dot.set_data_3d([spacecraft_pos[0]], [spacecraft_pos[1]], [spacecraft_pos[2]])
        spacecraft_visible = True
        
        # Update first transfer trail
        for i, pos in enumerate(spacecraft_pos):
            transfer_trail1_data[i].append(pos)
        transfer_trail1.set_data_3d(transfer_trail1_data[0], transfer_trail1_data[1], transfer_trail1_data[2])
    
    # Second transfer: Mars to Jupiter
    elif transfer2_start_day <= frame < transfer2_start_day + transfer2_duration_days:
        transfer_frame = frame - transfer2_start_day
        traj_idx = int((transfer_frame / transfer2_duration_days) * len(x_traj2))
        traj_idx = min(traj_idx, len(x_traj2) - 1)
        
        spacecraft_pos = [x_traj2[traj_idx], y_traj2[traj_idx], z_traj2[traj_idx]]
        spacecraft_dot.set_data_3d([spacecraft_pos[0]], [spacecraft_pos[1]], [spacecraft_pos[2]])
        spacecraft_visible = True
        
        # Update second transfer trail
        for i, pos in enumerate(spacecraft_pos):
            transfer_trail2_data[i].append(pos)
        transfer_trail2.set_data_3d(transfer_trail2_data[0], transfer_trail2_data[1], transfer_trail2_data[2])
    
    # Hide spacecraft when not in transfer
    if not spacecraft_visible:
        spacecraft_dot.set_data_3d([], [], [])
    
    return (earth_dot, mars_dot, jupiter_dot, spacecraft_dot, 
            earth_trail, mars_trail, jupiter_trail, 
            transfer_trail1, transfer_trail2)

# Create animation
total_frames = 4000

#Updates only every N frames
N = 5
anim = animation.FuncAnimation(fig, animate, frames=range(0, total_frames, N), interval=1, blit=False, repeat=False)

# Add legend and show
ax.legend()
plt.tight_layout()
plt.show()