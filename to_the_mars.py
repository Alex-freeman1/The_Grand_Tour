import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import lambertsolve_master as lt
import planetary_data as pd
from scipy.integrate import solve_ivp

n_steps = 4


spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")


def seconds(days):
    return days * 24 * 3600


def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

def norm(vec):
    return np.linalg.norm(vec)


def calc_vinfinity(tof, target, departure_planet, et0, vinf_target):
    """
    Function to compute difference between target vinfinity magnitude and
    current vinfinity magnitude for a given time of flight (tof).
    """
    # Position of planet1 at arrival time (et0 + tof)
    r1_planet1 = calc_ephemeris(target, et0 + tof, FRAME, OBSERVER)[:3]
    state_ofplanet0 = calc_ephemeris(departure_planet, et0, FRAME, OBSERVER)

    # Lambert solver: get spacecraft departure and arrival velocities    
    v0_sc_depart, v1_sc_arrive = lt.lambert_solver(
        state_ofplanet0[:3],    
        r1_planet1,   
        tof,                 # Time of flight in seconds
        sun_mu,               # Sun's gravitational parameter
        trajectory='pro'      # Prograde trajectory change to retro if doesn't work
    )

    # Compute vinfinity magnitude at departure planet
    vinf = np.linalg.norm(v0_sc_depart - state_ofplanet0[3:])

    # Return difference between target vinfinity and computed vinfinity
    return vinf_target - vinf

    
def vinfinity_match(planet0, planet1, v0_sc, et0, tof0,
                    diff_step=1e-3, tol=1e-4):
    """
    Given an incoming v-infinity vector to planet0, calculate the
    outgoing v-infinity vector that will arrive at planet1 after
    time of flight (tof) such that the incoming and outgoing v-infinity
    vectors at planet0 have equal magnitude.
    """
    
    # Get planet0 state at et0
    state0_planet0 = calc_ephemeris(planet0, et0, FRAME, OBSERVER)
    

    # Calculate magnitude of incoming V-infinity at planet0
    vinf = np.linalg.norm(v0_sc - state0_planet0[3:])

    # Define a local function wrapping calc_vinfinity for root finding
    def root_func(tof):
        return calc_vinfinity(tof, planet1, planet0, et0, vinf)

    # Solve for tof using Newton's method 
    tof, steps = lt.newton_root_single_fd(root_func, tof0, tol=tol, diff_step=diff_step)

    # Get position of planet1 at arrival
    r1_planet1 = calc_ephemeris(planet1, et0 + tof, FRAME, OBSERVER)

    # Compute Lambert solution for spacecraft velocities on transfer leg
    v_sc_depart, v_sc_arrive = lt.lambert_solver(
        state0_planet0[:3],    
        r1_planet1,   
        tof,                 # Time of flight in seconds
        sun_mu,               # Sun's gravitational parameter
        trajectory='pro'      # Prograde trajectory change to retro if doesn't work
    )
    return tof, v_sc_depart, v_sc_arrive


OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']      
                
# EVME
time0 = "1966-02-10"
time1 = "1966-07-07"
time2 = "1967-01-10"
time3= "1967-12-18"

planet0 = 'Earth' #Case sensitive
planet1 = 'Venus'
planet2 = 'Mars'

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
Venus = pd_req1['spice_name']
Mars = pd_req2['spice_name']

et_0 = spice.utc2et(time0)
et_1 = spice.utc2et(time1)

state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
state1 = calc_ephemeris(Venus, et_1, FRAME, OBSERVER)
tof = et_1 - et_0

vn_sc_depart, vn_sc_arrive = lt.lambert_solver(
    state0[:3],    
    state1[:3],   
    tof,                 # Time of flight in seconds
    sun_mu,               # Sun's gravitational parameter
    trajectory='pro'      # Prograde trajectory change to retro if doesn't work
)

periapsis = 0
turn_angle = 0
state_depart = np.concatenate((state0[:3], vn_sc_depart))
v_inf_0 = norm(vn_sc_depart - state0[3:])

state_sc_arrive = np.concatenate((state1[:3], vn_sc_arrive))
v_inf_1 = norm(vn_sc_arrive - state1[3:])


tof_guess = et_1 - et_0


tof, v_sc_depart, v_sc_arrive = vinfinity_match(
Earth, Venus,
state_sc_arrive[3:],  # incoming v at flyby planet
et_0, tof_guess)

state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
state1 = calc_ephemeris(Venus, et_0 + tof, FRAME, OBSERVER)
vinf_i = state_sc_arrive[3:] - state0[3:]
vinf_o = v_sc_depart - state0[3:]
tof_days = seconds(tof)
delta_angle = lt.vecs2angle(vinf_i, vinf_o) / 2.0


 