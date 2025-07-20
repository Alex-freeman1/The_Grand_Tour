# -*- coding: utf-8 -*-
"""
Created on Mon Jun 30 21:23:22 2025

@author: alexa
"""

import numpy as np
import spiceypy as spice
import datetime
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
# Your other functions (norm, radians, calc_ephemeris, rodrigues_rotate, etc.) go here

# Parameters
mu_sun = pd.sun['mu']
mu_mars = 42828.3
r_mars = 3396.2  # km
flyby_altitude = 300  # km
r_p = r_mars + flyby_altitude


def calc_ephemeris(target, ets, frame, observer):
    return np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])

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


def rodrigues_rotate(v, k, theta):
    return (v * np.cos(theta) +
            np.cross(k, v) * np.sin(theta) +
            k * np.dot(k, v) * (1 - np.cos(theta)))


from datetime import datetime, timedelta

def generate_date_range(start_date_str, end_date_str, step_days=1):
    start_date = datetime.strptime(start_date_str, "%Y-%m-%d")
    end_date = datetime.strptime(end_date_str, "%Y-%m-%d")
    delta = timedelta(days=step_days)
    
    dates = []
    current_date = start_date
    while current_date <= end_date:
        dates.append(current_date.strftime("%Y-%m-%d"))
        current_date += delta
    return dates

# Example usage:
earth_departure_dates = generate_date_range("2019-03-01", "2019-06-01", step_days=50)
mars_arrival_dates = generate_date_range("2019-03-01", "2025-06-01", step_days=50)
jupiter_arrival_dates = generate_date_range("2023-01-01", "2033-07-01", step_days=150)


best_match = None
best_score = np.inf

for dep_str in earth_departure_dates:
    for mars_arr_str in mars_arrival_dates:
        for jup_arr_str in jupiter_arrival_dates:
            try:
                # Convert dates to ET
                et_dep = spice.utc2et(dep_str)
                et_mars_arr = spice.utc2et(mars_arr_str)
                et_jup_arr = spice.utc2et(jup_arr_str)

                tof_earth_mars = et_mars_arr - et_dep
                tof_mars_jup = et_jup_arr - et_mars_arr
                if tof_earth_mars <= 0 or tof_mars_jup <= 0:
                    continue  # skip invalid times

                # Get planet states
                r_earth = calc_ephemeris(Earth, et_dep, FRAME, OBSERVER)[:3]
                v_earth = calc_ephemeris(Earth, et_dep, FRAME, OBSERVER)[3:]

                r_mars = calc_ephemeris(Mars, et_mars_arr, FRAME, OBSERVER)[:3]
                v_mars = calc_ephemeris(Mars, et_mars_arr, FRAME, OBSERVER)[3:]

                r_jup = calc_ephemeris(Jupiter, et_jup_arr, FRAME, OBSERVER)[:3]
                v_jup = calc_ephemeris(Jupiter, et_jup_arr, FRAME, OBSERVER)[3:]

                # Mars->Jupiter Lambert
                v_mars_depart, v_jup_arrive = lt.lambert_solver(
                    r_mars, r_jup, tof_mars_jup, mu_sun, trajectory='pro'
                )

                # Compute Mars hyperbolic excess velocity leaving Mars (relative velocity)
                v_inf_plus = v_mars_depart - v_mars
                v_inf_mag = np.linalg.norm(v_inf_plus)

                # Flyby turn angle delta
                ecc = 1 + (r_p * v_inf_mag**2) / mu_mars
                delta = 2 * np.arcsin(1 / ecc)

                # Rotate v_inf_plus by -delta to get incoming v_inf_minus
                dummy = np.array([0, 0, 1]) if abs(np.dot(v_inf_plus, [0, 0, 1])) < 0.99 else np.array([1, 0, 0])
                h_vec = np.cross(v_inf_plus, dummy)
                h_hat = h_vec / np.linalg.norm(h_vec)
                v_inf_minus = rodrigues_rotate(v_inf_plus, h_hat, -delta)

                # Incoming velocity at Mars (heliocentric)
                v_arrive_mars = v_inf_minus + v_mars

                # Compute spacecraft arrival position relative to Mars (heliocentric frame)
                v_hat = v_inf_minus / np.linalg.norm(v_inf_minus)
                r_rel = -r_p * v_hat
                r_arrival_sc = r_mars + r_rel

                # Earth->Mars Lambert with spacecraft arrival position
                v_earth_depart, v_mars_arrive = lt.lambert_solver(
                    r_earth, r_arrival_sc, tof_earth_mars, mu_sun, trajectory='pro'
                )

                # Now compare v_mars_arrive (from Earth->Mars leg) with v_arrive_mars (from flyby calc)
                delta_v = np.linalg.norm(v_mars_arrive - v_arrive_mars)

                # Save if this is the best match so far
                if delta_v < best_score:
                    best_score = delta_v
                    best_match = {
                        'earth_departure': dep_str,
                        'mars_arrival': mars_arr_str,
                        'jupiter_arrival': jup_arr_str,
                        'delta_v_error': delta_v,
                        'v_earth_depart': v_earth_depart,
                        'v_mars_arrive': v_mars_arrive,
                        'v_arrive_mars': v_arrive_mars,
                        'v_mars_depart': v_mars_depart,
                        'v_jup_arrive': v_jup_arrive,
                    }  

            except Exception as e:
                print(f"Error at {dep_str}, {mars_arr_str}, {jup_arr_str}: {e}")
                continue

print("Best match found:")
print(best_match)
