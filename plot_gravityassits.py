# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 13:21:32 2025

@author: alexa
"""

import sys, os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
#Import Python files
import matplotlib.pyplot as plt
import spiceypy as spice
import numpy as np
from datetime import datetime, timedelta


#Import Lambert tools
import lambertsolve_master as lt
import planetary_data as pd
import EMJ_gravityassist as ga

from joblib import Parallel, delayed
from tqdm import tqdm

# -------------------------------------------------
# Constants and SPICE setup
# -------------------------------------------------
OBSERVER = 'SUN'
FRAME = 'ECLIPJ2000'
mu_sun = 1.32712440018e11  # km^3/s^2, example value; adjust if you have pd.sun['mu']

# -------------------------------------------------
# Function to calculate ephemeris
# -------------------------------------------------
def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

# -------------------------------------------------
# Johann function
# -------------------------------------------------
def johann(dep_planet, arr_planet, departure0, departure1, arrival0, arrival1, cutoff_c3=200):
    """Compute C3 shorts and spacecraft departure velocities for a porkchop grid"""
    
    # Step size
    step_size = 50 if dep_planet in ["EARTH", "MARS BARYCENTER"] else 60
    step = step_size * 24*3600  # seconds

    # Time arrays
    et_departures = np.arange(spice.utc2et(departure0), spice.utc2et(departure1)+step, step)
    et_arrivals   = np.arange(spice.utc2et(arrival0), spice.utc2et(arrival1)+step, step)
    ds = len(et_departures)
    as_ = len(et_arrivals)

    # Output arrays
    C3_shorts = np.zeros((as_, ds))
    vsc_departures = np.zeros((as_, ds, 3))
    vsc_arrivals = np.zeros((as_, ds, 3))
    
    # Precompute ephemerides
    ephem_departures = [calc_ephemeris(dep_planet, et, FRAME, OBSERVER) for et in et_departures]
    ephem_arrivals   = [calc_ephemeris(arr_planet, et, FRAME, OBSERVER) for et in et_arrivals]

    # Lambert solver for each combination
    def compute_lambert_entry(states_depart, states_arrive, tof, sun_mu, cutoff_c3):
        from numpy.linalg import norm
        if tof <= 0 or norm(states_depart[:3] - states_arrive[:3]) < 1e6:
            return cutoff_c3, np.zeros(3), np.zeros(3)
        try:
            v_sc_depart, v_sc_arrive = lt.lambert_solver(states_depart[:3], states_arrive[:3], tof, sun_mu, trajectory='pro')
        except:
            v_sc_depart, v_sc_arrive = states_depart[3:] + 1000, states_depart[3:] + 1000
        C3_short = min(norm(v_sc_depart - states_depart[3:])**2, cutoff_c3)
        return C3_short, v_sc_depart, v_sc_arrive

    # Parallel computation
    results = Parallel(n_jobs=-1)(
        delayed(compute_lambert_entry)(
            ephem_departures[nd],
            ephem_arrivals[na],
            et_arrivals[na]-et_departures[nd],
            mu_sun,
            cutoff_c3
        )
        for na in tqdm(range(as_), desc="Calculating Transfer")
        for nd in range(ds)
    )

    # Fill output arrays
    for idx, (C3_short, v_sc_depart, v_sc_arrive) in enumerate(results):
        na = idx // ds
        nd = idx % ds
        C3_shorts[na, nd] = C3_short
        vsc_departures[na, nd] = v_sc_depart
        vsc_arrivals[na, nd] = v_sc_arrive

    # Normalize departure/arrival times in days
    normed_departures = (et_departures - et_departures[0]) / (3600*24)
    normed_arrivals   = (et_arrivals - et_arrivals[0]) / (3600*24)

    return normed_departures, normed_arrivals, C3_shorts, vsc_departures, vsc_arrivals

# -------------------------------------------------
# V-infinity matching for gravity assists
# -------------------------------------------------
def calc_vinfinity(tof, target_planet, departure_planet, et0, vinf_target):
    r1_target = calc_ephemeris(target_planet, et0 + tof, FRAME, OBSERVER)[:3]
    state_departure = calc_ephemeris(departure_planet, et0, FRAME, OBSERVER)
    v0_sc_depart, _ = lt.lambert_solver(state_departure[:3], r1_target, tof, mu_sun, trajectory='pro')
    vinf = np.linalg.norm(v0_sc_depart - state_departure[3:])
    return vinf_target - vinf

def vinfinity_match(planet0, planet1, v_sc_incoming, et0, tof0, diff_step=1e-3, tol=1e-4):
    state0_planet0 = calc_ephemeris(planet0, et0, FRAME, OBSERVER)
    vinf_incoming = v_sc_incoming - state0_planet0[3:]
    vinf_magnitude = np.linalg.norm(vinf_incoming)

    def root_func(tof):
        return calc_vinfinity(tof, planet1, planet0, et0, vinf_magnitude)

    tof, _ = lt.newton_root_single_fd(root_func, tof0, tol=tol, diff_step=diff_step)
    
    r1_planet1 = calc_ephemeris(planet1, et0 + tof, FRAME, OBSERVER)[:3]
    v_sc_depart, v_sc_arrive = lt.lambert_solver(state0_planet0[:3], r1_planet1, tof, mu_sun, trajectory='pro')

    return tof, v_sc_depart, v_sc_arrive

# -------------------------------------------------
# Earth -> Jupiter -> Saturn
# -------------------------------------------------

pd_bodies = pd.bodies
def loop_bodies(planeti, planetf):
    pd_req0 = pd_req1 = None
    for body in pd_bodies:
        if body['name'] == planeti:
            pd_req0 = body['spice_name']
        elif body['name'] == planetf:
            pd_req1 = body['spice_name']
        
    return pd_req0, pd_req1


# Set up planets
planet0 = 'Earth'
planet1 = 'Jupiter'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)


time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival


# Earth -> Jupiter porkchop
dep_dates_0 = '1976-02-01'
dep_dates_1 = '1979-11-07'
arr_dates_0 = '1978-02-01'
arr_dates_1 = '1980-12-01'

norm_dep, norm_arr, C3_shorts, departure_velocities, arrival_velocities = johann(
    departure_planet, arrival_planet,
    dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1
)

# Jupiter -> Saturn
planet0 = 'Jupiter'
planet1 = 'Saturn'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_2_i = "1978-02-01"
dep_dates_2_f = "1980-12-01"
arr_dates_2_i = "1980-01-01"
arr_dates_2_f = "1983-01-01"

# Compute ET arrays for second leg
step = 50*24*3600  # example


et_departures = np.arange(spice.utc2et(dep_dates_2_i), spice.utc2et(dep_dates_2_f)+step, step)
et_arrivals   = np.arange(spice.utc2et(arr_dates_2_i), spice.utc2et(arr_dates_2_f)+step, step)

print(et_departures)
ds = len(et_departures)
as_ = len(et_arrivals)

# Loop over all departure/arrival combinations

q = 1


success_dep = []
success_arr = []
for na in range(as_):
    for nd in range(ds):
        tof = et_arrivals[na] - et_departures[nd]
        prev_na = min(na, arrival_velocities.shape[0]-1)
        prev_nd = min(nd, arrival_velocities.shape[1]-1)
        v_sc_incoming = arrival_velocities[prev_na, prev_nd]

        try:
            tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
                departure_planet, arrival_planet,
                v_sc_incoming,
                et_departures[nd], tof
            )
            #print(f"Success: tof_2 = {tof_2} seconds for departure index {nd}, arrival index {na}")
            # Departure date
            print(q)
            dep_date_utc = spice.et2utc(et_departures[nd], 'C', 0)
            print(f"Departure: {dep_date_utc}")
            # Arrival date = original arrival + new TOF
            arr_date_utc = spice.et2utc(et_departures[nd] + tof_2, 'C', 0)
            print(f"Arrival (with new TOF): {arr_date_utc}")
            
            dep_date = datetime.strptime(dep_date_utc, "%Y %b %d %H:%M:%S")
            arr_date = datetime.strptime(arr_date_utc, "%Y %b %d %H:%M:%S")
            
            success_dep.append(dep_date)
            success_arr.append(arr_date)


            q = q + 1
        except RuntimeError as e:
            #print(f"Warning: Newton iteration failed for departure index {nd}, arrival index {na}. Skipping.")
            continue  # skip this case
            


# Define linewdith
lw = 1.5

'''
Create The Plots
'''
fig, ax = plt.subplots()
ax.scatter(success_dep, success_arr, color='m', s=50)  # s = marker size
ax.set_xlabel("Departure Date")
ax.set_ylabel("Arrival Date")
ax.set_title("Successful Lambert Solution")

# Improve date formatting
fig.autofmt_xdate()
plt.show()