
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
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
# time0 = "1977-08-15"  # Earth departure
# time1 = "1979-07-22"  # Jupiter flyby

time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival


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
print(f"Heliocentric arriving veloicty: {v_sc_arrive_1} km/s)")
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
    
     





     
