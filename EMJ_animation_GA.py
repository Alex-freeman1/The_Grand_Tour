# -*- coding: utf-8 -*-
"""
Created on Mon Jul  7 17:20:02 2025

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

def seconds(days):
    return days * 24 * 3600

def calc_ephemeris(target, ets, frame, observer):
    return np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])


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



# Dates for Planetary flybys

time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival


start = "1977-02-15"

et_0 = spice.utc2et(time0)
et_1 = spice.utc2et(time1)
et_2 = spice.utc2et(time2)


# Earth to Jupiter time
tof_1 = et_1 - et_0   
      

print(f"Earth to Jupiter: {tof_1/86400:.1f} days")
#print(f"Jupiter to Saturn: {tof2/86400:.1f} days")

# Get ephemeris data for all mission points
state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
state1 = calc_ephemeris(Jupiter, et_1, FRAME, OBSERVER)

# Solve Lambert's problem for Earth to Jupiter
try:
    v_sc_depart_1, v_sc_arrive_1 = lt.lambert_solver(
        state0[:3],    
        state1[:3],   
        tof_1,
        mu_sun,
        trajectory='pro'
    )

except Exception as e:
    print(f"Earth-Mars Lambert solver failed: {e}")
    v_sc_depart_1 = np.array([1000, 1000, 1000])
    v_sc_depart_1 = np.array([1000, 1000, 1000])



v_inf_dep_E = v_sc_depart_1 - state0[3:]
v_inf_arr_J = v_sc_arrive_1 - state1[3:]

tof_guess_2 = et_2 - et_1

"""
This function takes the input as the heliocentric arrival velocity to Jupiter (from Earth)
and using the tof_guess trys to find the heliocentric leaving velocity such that it will both get to 
Saturn and that the v_inf both plus and minus are both equal in magnitude. Once this is solved, it returns 
the new time of flight as well as the departing heliocentric velocity
"""
tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
    Jupiter, Saturn,
    v_sc_arrive_1,  # spacecraft velocity at Jupiter arrival
    et_1, tof_guess_2
)



print(v_sc_depart_2)

vinf_in = v_sc_arrive_1 - state1[3:]
vinf_out = v_sc_depart_2 - state1[3:]
def_angle = np.arccos(np.dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))
r_jup = 69911
mu_jup = 126.687*10**6


def altitude(delta_a, v_infx):
    e = 1 / np.sin(delta_a / 2)
    r_p = (e - 1) * mu_jup / v_infx**2
    h = r_p - r_jup
    return h

v_inf_mag = norm(vinf_in)
h = altitude(def_angle, v_inf_mag)
print(h)

X0 = [*state0[:3], *v_sc_depart_1]      # Earth to Jupiter
X1 = [*state1[:3], *v_sc_depart_2]       # Jupiter to Saturn

def two_body(t, y):
    """Two-body dynamics for spacecraft motion"""
    r = y[:3]
    v = y[3:]
    norm_r = np.linalg.norm(r)
    a = -mu_sun * r / norm_r**3
    return np.concatenate((v, a))



t_span1 = (0, tof_1)
t_eval1 = np.linspace(0, tof_1, 1000)
solution1 = solve_ivp(two_body, t_span1, X0, t_eval=t_eval1, method='RK45', rtol=1e-12, atol=1e-12)

t_span2 = (0, tof_2)
t_eval2 = np.linspace(0, tof_2, 1000)
solution2 = solve_ivp(two_body, t_span2, X1, t_eval=t_eval2, method='RK45', rtol=1e-12, atol=1e-12)

# Extract trajectory data
if solution1.success:
    x_traj1, y_traj1, z_traj1 = solution1.y[0], solution1.y[1], solution1.y[2]
else:
    x_traj1 = y_traj1 = z_traj1 = np.array([0])
    print("First trajectory integration failed")

if solution2.success:
    x_traj2, y_traj2, z_traj2 = solution2.y[0], solution2.y[1], solution2.y[2]
else:
    x_traj2 = y_traj2 = z_traj2 = np.array([0])
    print("Second trajectory integration failed")


et_start = int(spice.utc2et(start))
time_step = 86400  # 1 day


def time_orbit(a, μ):
    T = 2 * np.pi * np.sqrt(a**3 / μ)
    time = T 
    return time


# Semi-major axes (km)
a_earth   = 149597870.7
a_jupiter = 778547200.0
a_saturn  = 1433449370.0

T_earth   = time_orbit(a_earth, mu_sun) / 86400  # in days
T_jupiter = time_orbit(a_jupiter, mu_sun) / 86400
T_saturn  = time_orbit(a_saturn, mu_sun) / 86400

time_step = 1 * 86400  # 1 day in seconds

def planet_positions(spk_id, et_start, T_days, step_days=1):
    times = np.arange(et_start, et_start + T_days * 86400, step_days * 86400)
    return [calc_ephemeris(spk_id, t, FRAME, OBSERVER)[:3] for t in times]

planet_0   = planet_positions(Earth,   et_start, T_earth)
planet_1 = planet_positions(Jupiter, et_start, T_jupiter)
planet_2  = planet_positions(Saturn,  et_start, T_saturn)


# FIXED: Calculate transfer timing in animation frames
transfer1_start_frame = int((et_0 - et_start) / time_step)
transfer1_duration_frames = int(tof_1 / time_step)
transfer2_start_frame = int((et_1 - et_start) / time_step)
transfer2_duration_frames = int(tof_2 / time_step)

print(f"Transfer 1: frames {transfer1_start_frame} to {transfer1_start_frame + transfer1_duration_frames}")
print(f"Transfer 2: frames {transfer2_start_frame} to {transfer2_start_frame + transfer2_duration_frames}")

# Set up the 3D figure and axis
fig = plt.figure(figsize=(14, 10))
ax = fig.add_subplot(111, projection='3d')

# Set axis limits
margin = 1.2
max_range = max(np.max(np.abs(planet_2)), np.max(np.abs(planet_1)), np.max(np.abs(planet_0)))
ax.set_xlim(-max_range * margin, max_range * margin)
ax.set_ylim(-max_range * margin, max_range * margin)
ax.set_zlim(-max_range * margin, max_range * margin)

# Labels and title
ax.set_xlabel('X (km)')
ax.set_ylabel('Y (km)')
ax.set_zlabel('Z (km)')
ax.set_title('Earth-Mars-Jupiter Transfer Trajectory')

# Create plot objects
zero_dot, = ax.plot([], [], [], 'bo', markersize=8, label='Earth')
one_dot, = ax.plot([], [], [], 'ro', markersize=8, label='Jupiter')
two_dot, = ax.plot([], [], [], 'orange', marker='o', markersize=10, label='Saturn')
spacecraft_dot, = ax.plot([], [], [], 'go', markersize=6, label='Spacecraft')
sun_dot, = ax.plot([0], [0], [0], 'yo', markersize=12, label='Sun')

# Trail objects
zero_trail, = ax.plot([], [], [], 'b-', alpha=0.3, linewidth=1)
one_trail, = ax.plot([], [], [], 'r-', alpha=0.3, linewidth=1)
two_trail, = ax.plot([], [], [], 'orange', alpha=0.3, linewidth=1)
transfer_trail1, = ax.plot([], [], [], 'g-', alpha=0.8, linewidth=2, label='Earth-Jupiter')
transfer_trail2, = ax.plot([], [], [], 'm-', alpha=0.8, linewidth=2, label='Jupiter-Saturn')

# Trail data
trail_length = 50
zero_trail_data = [[], [], []]
one_trail_data = [[], [], []]
two_trail_data = [[], [], []]
transfer_trail1_data = [[], [], []]
transfer_trail2_data = [[], [], []]

def animate(frame):
    """Animation function"""
    # Ensure frame is within bounds
    #frame = frame % len(earth_pos)
    
    frame_earth   = frame % len(planet_0)
    frame_mars    = frame % len(planet_1)
    frame_jupiter = frame % len(planet_2)
    
    # Current planet positions
    zero_current = planet_0[frame_earth]
    one_current = planet_1[frame_mars]
    two_current = planet_2[frame_jupiter]
    
    # Update planet positions
    zero_dot.set_data_3d([zero_current[0]], [zero_current[1]], [zero_current[2]])
    one_dot.set_data_3d([one_current[0]], [one_current[1]], [one_current[2]])
    two_dot.set_data_3d([two_current[0]], [two_current[1]], [two_current[2]])
    
    # Update trails every few frames
    if frame % 3 == 0:
        for i, pos in enumerate(zero_current):
            zero_trail_data[i].append(pos)
            if len(zero_trail_data[i]) > trail_length:
                zero_trail_data[i].pop(0)
        
        for i, pos in enumerate(one_current):
            one_trail_data[i].append(pos)
            if len(one_trail_data[i]) > trail_length:
                one_trail_data[i].pop(0)
        
        for i, pos in enumerate(two_current):
            two_trail_data[i].append(pos)
            if len(two_trail_data[i]) > trail_length:
                two_trail_data[i].pop(0)
    
    zero_trail.set_data_3d(zero_trail_data[0], zero_trail_data[1], zero_trail_data[2])
    one_trail.set_data_3d(one_trail_data[0], one_trail_data[1], one_trail_data[2])
    two_trail.set_data_3d(two_trail_data[0], two_trail_data[1], two_trail_data[2])
    
    
    spacecraft_visible = False
    
    # First transfer timing
    if transfer1_start_frame <= frame < transfer1_start_frame + transfer1_duration_frames:
        transfer_progress = (frame - transfer1_start_frame) / transfer1_duration_frames
        traj_idx = int(transfer_progress * (len(x_traj1) - 1))
        traj_idx = min(traj_idx, len(x_traj1) - 1)
        
        spacecraft_pos = [x_traj1[traj_idx], y_traj1[traj_idx], z_traj1[traj_idx]]
        spacecraft_dot.set_data_3d([spacecraft_pos[0]], [spacecraft_pos[1]], [spacecraft_pos[2]])
        spacecraft_visible = True
        
        # Update first transfer trail
        for i, pos in enumerate(spacecraft_pos):
            transfer_trail1_data[i].append(pos)
        transfer_trail1.set_data_3d(transfer_trail1_data[0], transfer_trail1_data[1], transfer_trail1_data[2])
    
    # Second transfer timing
    elif transfer2_start_frame <= frame < transfer2_start_frame + transfer2_duration_frames:
        transfer_progress = (frame - transfer2_start_frame) / transfer2_duration_frames
        traj_idx = int(transfer_progress * (len(x_traj2) - 1))
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
    
    return (zero_dot, one_dot, two_dot, spacecraft_dot, 
            zero_trail, one_trail, two_trail, 
            transfer_trail1, transfer_trail2)

# Create animation

total_frames = transfer2_start_frame + transfer2_duration_frames + 200  

# Update every N frames to speed up animation
N = 5
anim = animation.FuncAnimation(fig, animate, frames=range(0, total_frames, N), 
                              interval=10, blit=False, repeat=False)

# Add legend and show
ax.legend()
plt.tight_layout()
plt.show()