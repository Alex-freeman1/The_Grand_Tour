# -*- coding: utf-8 -*-
"""
Created on Thu May 29 23:06:48 2025

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

pd_bodies = pd.bodies

for i in range(10):
    if pd_bodies[i]['name'] == planet0:
        pd_req0 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet1:
        pd_req1 = pd_bodies[i]

# Holds the spice name of the planets        
Earth = pd_req0['spice_name']
Mars = pd_req1['spice_name']


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
departure0 = '2006-02-15'
arrival0 = '2007-11-10'  

et_departure = spice.utc2et(departure0)
et_arrival = spice.utc2et(arrival0)

tof1 = et_arrival - et_departure

# Get ephemeris data for departure and arrival
ephem_departure = calc_ephemeris(Earth, et_departure, FRAME, OBSERVER)
ephem_arrival = calc_ephemeris(Mars, et_arrival, FRAME, OBSERVER)




# Solve Lambert's problem
try:
    v_sc_depart, v_sc_arrive = lt.lambert_solver(
        ephem_departure[:3],  # Earth position at departure
        ephem_arrival[:3],    # Mars position at arrival
        tof1,                 # Time of flight in seconds
        mu_sun,               # Sun's gravitational parameter
        trajectory='pro'      # Prograde trajectory
    )

except Exception as e:
    print(f"Lambert solver failed: {e}")
    v_sc_depart = np.array([1000, 1000, 1000])
    v_sc_arrive = np.array([1000, 1000, 1000])


vx0, vy0, vz0 = v_sc_depart


# Initial conditions for spacecraft trajectory
x0, y0, z0 = ephem_departure[:3] 


X0 = [x0, y0, z0, vx0, vy0, vz0] 


def two_body(t, y):
    """Two-body dynamics for spacecraft motion"""
    r = y[:3]
    v = y[3:]
    norm_r = np.linalg.norm(r)
    a = -mu_sun * r / norm_r**3
    return np.concatenate((v, a))

# Integrate spacecraft trajectory
t_span = (0, tof1)
t_eval = np.linspace(0, tof1, 1000)

solution = solve_ivp(
    two_body, t_span, X0, t_eval=t_eval, method='RK45',
    rtol=1e-12, atol=1e-12)



if solution.success:
    x_traj, y_traj, z_traj = solution.y[0], solution.y[1], solution.y[2]
        
else:
    x_traj = y_traj = z_traj = np.array([0])
    



# Setup orbital data for planets
start = '2005-11-26' 
et_start = int(spice.utc2et(start))
earth_orbital_period = seconds(365.25)
mars_orbital_period = seconds(687)


time_step = 86400  # 1 day

# Calculate when transfer begins and ends in simulation time
transfer_start_day = int((et_departure - et_start) / time_step)
transfer_duration_days = int(tof1 / time_step)


# Generate planet positions
mars_pos = []
earth_pos = []

et_end = et_start + int(earth_orbital_period) 
et_end_mars = et_start + int(mars_orbital_period) 

for t in range(et_start, et_end, time_step):  # Use time_step instead of every second
    earth_pos.append(calc_ephemeris(Earth, t, FRAME, OBSERVER)[:3])
    
for t1 in range(et_start, et_end_mars, time_step):
    mars_pos.append(calc_ephemeris(Mars, t1, FRAME, OBSERVER)[:3])

mars_pos = np.array(mars_pos)
earth_pos = np.array(earth_pos)

# Set up the 3D figure and axis
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')


margin = 1.2
max_range = max(np.max(np.abs(mars_pos)), np.max(np.abs(earth_pos)))
ax.set_xlim(-max_range * margin, max_range * margin)
ax.set_ylim(-max_range * margin, max_range * margin)
ax.set_zlim(-max_range * margin, max_range * margin)

# Labels and title
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('Earth-Mars Transfer Trajectory')

# Create plot objects
earth_dot, = ax.plot([], [], [], 'bo', markersize=8, label='Earth')
mars_dot, = ax.plot([], [], [], 'ro', markersize=8, label='Mars')
spacecraft_dot, = ax.plot([], [], [], 'go', markersize=6, label='Spacecraft')
sun_dot, = ax.plot([0], [0], [0], 'yo', markersize=12, label='Sun')

# Trail objects
earth_trail, = ax.plot([], [], [], 'b-', alpha=0.3, linewidth=1)
mars_trail, = ax.plot([], [], [], 'r-', alpha=0.3, linewidth=1)
transfer_trail, = ax.plot([], [], [], 'g-', alpha=0.8, linewidth=2)

# Trail data
trail_length = 1000
earth_trail_data = [[], [], []]
mars_trail_data = [[], [], []]
transfer_trail_data = [[], [], []]

def animate(frame):
    """Animation function"""
    # Handle wraparound for long simulations
    earth_idx = frame % len(earth_pos)
    mars_idx = frame % len(mars_pos)
    
    # Current planet positions
    earth_current = earth_pos[earth_idx]
    mars_current = mars_pos[mars_idx]
    
    # Update planet positions
    earth_dot.set_data_3d([earth_current[0]], [earth_current[1]], [earth_current[2]])
    mars_dot.set_data_3d([mars_current[0]], [mars_current[1]], [mars_current[2]])
    
    # Update planet trails
    for i, pos in enumerate(earth_current):
        earth_trail_data[i].append(pos)
        if len(earth_trail_data[i]) > trail_length:
            earth_trail_data[i].pop(0)
    
    for i, pos in enumerate(mars_current):
        mars_trail_data[i].append(pos)
        if len(mars_trail_data[i]) > trail_length:
            mars_trail_data[i].pop(0)
    
    earth_trail.set_data_3d(earth_trail_data[0], earth_trail_data[1], earth_trail_data[2])
    mars_trail.set_data_3d(mars_trail_data[0], mars_trail_data[1], mars_trail_data[2])
    
    # Handle spacecraft trajectory during transfer window
    spacecraft_visible = False
    if transfer_start_day <= frame < transfer_start_day + transfer_duration_days:
        # We're in the transfer window
        transfer_frame = frame - transfer_start_day
        traj_idx = int((transfer_frame / transfer_duration_days) * len(x_traj))
        traj_idx = min(traj_idx, len(x_traj) - 1)
        
        spacecraft_pos = [x_traj[traj_idx], y_traj[traj_idx], z_traj[traj_idx]]
        spacecraft_dot.set_data_3d([spacecraft_pos[0]], [spacecraft_pos[1]], [spacecraft_pos[2]])
        spacecraft_visible = True
        
        # Update transfer trail
        for i, pos in enumerate(spacecraft_pos):
            transfer_trail_data[i].append(pos)
        
        transfer_trail.set_data_3d(transfer_trail_data[0], transfer_trail_data[1], transfer_trail_data[2])
    
    # Hide spacecraft when not in transfer
    if not spacecraft_visible:
        spacecraft_dot.set_data_3d([], [], [])
    
    return earth_dot, mars_dot, spacecraft_dot, earth_trail, mars_trail, transfer_trail

# Create animation
total_frames = 2000
anim = animation.FuncAnimation(fig, animate, frames=total_frames, 
                             interval=50, blit=False, repeat=False)

# Add legend and show
ax.legend()
plt.tight_layout()
plt.show()

# Optional: Save animation (uncomment if needed)
# anim.save('earth_mars_transfer.gif', writer='pillow', fps=20)