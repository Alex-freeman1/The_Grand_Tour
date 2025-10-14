# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 13:21:32 2025

@author: alexa
"""

import sys, os

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


# Establish Kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")




OBSERVER = pd.sun['name']
FRAME = 'ECLIPJ2000'
mu_sun = pd.sun['mu']

def norm(vec):
    return np.linalg.norm(vec)

def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])


def perihelion_distance(r1, v1, mu):
    
    h_vec = np.cross(r1, v1)
    h = norm(h_vec)
    energy = norm(v1)**2 / 2 - mu / norm(r1)
    a = -mu / (2 * energy)
    e = np.sqrt(1 - h**2 / (a * mu))
    rp = a * (1 - e)
    return rp



def johann(dep_planet, arr_planet, departure0, departure1, arrival0, arrival1, cutoff_c3=200):
    
    # Step size
    step_size = 50
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
            
            rp = perihelion_distance(states_depart[:3], v_sc_depart, sun_mu)
            if rp < 0.8 * 1.496e8:  
                return cutoff_c3, np.zeros(3), np.zeros(3)
        
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

# Earth -> Jupiter porkchop
dep_dates_0 = '1977-02-01'
dep_dates_1 = '1978-11-07'
arr_dates_0 = '1978-02-01'
arr_dates_1 = '1980-12-01'


planet_0_dep = spice.utc2et(dep_dates_0)
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
step = 50*24*3600  

et_departures = np.arange(spice.utc2et(dep_dates_2_i), spice.utc2et(dep_dates_2_f)+step, step)
et_arrivals   = np.arange(spice.utc2et(arr_dates_2_i), spice.utc2et(arr_dates_2_f)+step, step)



#print(et_departures)
ds = len(et_departures)
as_ = len(et_arrivals)


r_jup = 69911
mu_jup = 126.687*10**6

def altitude(delta_a, v_infx):
    e = 1 / np.sin(delta_a / 2)
    r_p = (e - 1) * mu_jup / v_infx**2
    h = r_p - r_jup
    return h

terra_dates = []
jup_dates = []
sat_dates = []

for i in range(arrival_velocities.shape[1]): 

    for na in range(as_):
        for nd in range(ds):
            
            
            tof = et_arrivals[na] - et_departures[nd]

            v_sc_incoming = arrival_velocities[nd, i]

            try:
                tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
                    departure_planet, arrival_planet,
                    v_sc_incoming,
                    et_departures[nd], tof
                )
                
            
            
            except RuntimeError:
                #print(f"Warning: Newton iteration failed for departure index {nd}, arrival index {na}. Skipping.")
                continue  # skip this case
                
            else:   
                
                state_Jupiter_flyby = calc_ephemeris(departure_planet, et_departures[nd], FRAME, OBSERVER)
                dep_ter = planet_0_dep + (i * 50 * 86400)
                dep_jup = et_departures[nd]
                arr_sat = et_departures[nd] + tof_2
                
                
                vinf_in = v_sc_incoming - state_Jupiter_flyby[3:]
                vinf_out = v_sc_depart_2 - state_Jupiter_flyby[3:]
                    
                def_angle = np.arccos(np.dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))
                   
                v_inf_mag = norm(vinf_in)
                h = altitude(def_angle, v_inf_mag)
                
                if h > 0 and vinf_in[0] < 0:
                    
                    terra_dates.append(dep_ter)
                    jup_dates.append(dep_jup)
                    sat_dates.append(arr_sat)
                    
                    print(" --- Success --- ")
                    # print(nd, na)
                    print(f"et_0 = {dep_ter}")
                    print(f"et_1 = {dep_jup}")
                    print(f"et_2 = {arr_sat}")
                    # print('\n')
                    # print(f"The altitude is {h} km")
                    # print(f"incoming: {v_sc_incoming}")
                    # print(f"Outgoing: {v_sc_depart_2}")
                    # print(vinf_in)
                    # print(vinf_out)
                    # print('\n')
                    # print('\n')
            


# # Define linewdith
# lw = 1.5

# '''
# Create The Plots
# '''
# fig, ax = plt.subplots()
# ax.scatter(success_dep, success_arr, color='m', s=50)  # s = marker size
# ax.set_xlabel("Departure Date")
# ax.set_ylabel("Arrival Date")
# ax.set_title("Successful Lambert Solution")

# # Improve date formatting
# fig.autofmt_xdate()
# plt.show()