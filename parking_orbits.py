# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 20:42:26 2024

@author: alexa
"""
"""
If error "Cannot find spicepy or load kernels"
Ensure the script is run in all folders (top right) so that the file can find relevant data
"""
#Import Python files
import matplotlib.pyplot as plt
import spiceypy as spice

import numpy as np

#Import Lambert tools
import lambertSolver as lt
import planetary_data as pd

from joblib import Parallel, delayed
from tqdm import tqdm

# plt.style.use( 'dark_background' )



'''
Setting up Ephemeris Data and Kernels
'''



# Establish Kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")


# Input code - this is variable, please change to whatever you want

planet0 = 'Earth' #Case sensitive
planet1 = 'Jupiter'

#year / month / day
departure0 = '2006-02-01'         # Intial departure date
arrival0 = '2009-01-01'        # Initial arrival date
         
# --------------------



# Frame of reference from the planetary data
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']





# Loops through the podies to check against the inputted planet0 and planet1
pd_bodies = pd.bodies
for i in range(10):
    if pd_bodies[i]['name'] == planet0:
        pd_req0 = pd_bodies[i]
    elif pd_bodies[i]['name'] == planet1:
        pd_req1 = pd_bodies[i]
        
# Holds the spice name of the planets        
departure_planet = pd_req0['spice_name']
arrival_planet = pd_req1['spice_name']



#cut off velocity (if the lambert solver fails to converge)
cutoff_v = 100
cutoff_c3 = cutoff_v ** 2


#Fetch the time data and store into arrays from initial to final
et_departure = spice.utc2et(departure0)
et_arrival   = spice.utc2et(arrival0)
           
tof = et_arrival - et_departure


'''
Using Lambert's Solver to return the required velocities
'''

          

earth_rad = 6378.137  # km
h_alt = 1000

def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

def calc_ephemerisE(target, ets, frame, observer, theta):
    state = np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])
    pos = state[:3]
    vel = state[3:]
    
        
    r_offset = earth_rad + h_alt
        
    #xy Plane
    #offset = np.array([np.cos(theta), np.sin(theta), 0.0]) * r_offset
    
    #yz Plane
    # offset = np.array([0.0, np.cos(theta), np.sin(theta)])* r_offset
    
    #xz Plane
    offset = np.array([np.cos(theta), 0.0, np.sin(theta)])* r_offset
    
    r_sc = pos + offset

    return np.concatenate([r_sc, vel])

	
def norm(vec):
    return np.linalg.norm(vec)

n_angles = 1000
angles = np.linspace(0.0, 2*np.pi, n_angles, endpoint=False)


states_arrive = calc_ephemeris(arrival_planet, et_arrival, FRAME, OBSERVER)


energies_parking = []
for psi in angles:
    states_depart = calc_ephemerisE(departure_planet, et_departure, FRAME, OBSERVER, psi) 
    v_sc_depart_short, _ = lt.lambert_solver(states_depart[:3], states_arrive[:3], tof, sun_mu, trajectory='pro')
    C3_short = min(norm(v_sc_depart_short - states_depart[3:])**2, cutoff_c3)
    
    energies_parking.append(C3_short)
#C3_long = min(norm(v_sc_depart_long - states_depart[3:])**2, cutoff_c3)


plt.plot(angles, energies_parking)
plt.show()



