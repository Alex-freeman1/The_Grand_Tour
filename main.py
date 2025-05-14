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
planet1 = 'Mars'

departure0 = '2026-06-20'         # Intial departure date
departure1 = '2030-11-07'         # Final departure date
arrival0 = '2026-12-01'         # Initial arrival date
arrival1 = '2031-02-24'          # Final arrival date


# Step size (in days)
step_size = 10

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

#Convert from days into seconds
step = step_size*3600*24

#Fetch the time data and store into arrays from initial to final
et_departures = np.arange(spice.utc2et(departure0), spice.utc2et(departure1) + step, step)
et_arrivals   = np.arange(spice.utc2et(arrival0), spice.utc2et(arrival1) + step, step)
           
#Find the length of these time arrays
ds = len(et_departures)
as_ = len(et_arrivals)

total = as_ * ds   


#Set up zerod-arrays to store values in the future
C3_shorts     = np.zeros( (as_, ds) )
C3_longs      = np.zeros( (as_, ds) )
# v_inf_shorts  = np.zeros( (as_, ds) )
# v_inf_longs   = np.zeros( (as_, ds) )
tofs          = np.zeros( (as_, ds) )  

#Set up arrays of length equal to launchwindows 
x = np.arange(ds)
y = np.arange(as_)   
     

'''
Using Lambert's Solver to return the required velocities
'''

          

#Define function to get the emphemeris data - array holding positions and velocity vectors
def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

	
def norm(vec):
    return np.linalg.norm(vec)


# Retrieve the position and velocity vectors and store in a sixth-lengthed vector from function above
ephem_departures = [calc_ephemeris(departure_planet, et, FRAME, OBSERVER) for et in et_departures]
ephem_arrivals = [calc_ephemeris(arrival_planet, et, FRAME, OBSERVER) for et in et_arrivals]
          

    

def compute_lambert_entry(na, nd):
    #print( f'{na + 1} / {as_}.' )
    states_depart = ephem_departures[nd]
    states_arrive = ephem_arrivals[na]
    tof = et_arrivals[na] - et_departures[nd]

    if tof <= 0 or norm(states_depart[:3] - states_arrive[:3]) < 1e6:
        return cutoff_c3, cutoff_c3, tof

    try:
        v_sc_depart_short, _ = lt.lambert_solver(states_depart[:3], states_arrive[:3], tof, sun_mu, trajectory='pro')
    except:
        v_sc_depart_short = states_depart[3:] + 1000

    try:
        v_sc_depart_long, _ = lt.lambert_solver(states_depart[:3], states_arrive[:3], tof, sun_mu, trajectory='retro')
    except:
        v_sc_depart_long = states_depart[3:] + 1000

    C3_short = min(norm(v_sc_depart_short - states_depart[3:])**2, cutoff_c3)
    C3_long = min(norm(v_sc_depart_long - states_depart[3:])**2, cutoff_c3)
    
    return C3_short, C3_long, tof


results = Parallel(n_jobs=-1)(
    delayed(compute_lambert_entry)(na, nd)
    for na in tqdm(range(as_), desc="na loop")
    for nd in range(ds)
    
)

for idx, (C3_short, C3_long, tof) in enumerate(results):
    na = idx // ds
    nd = idx % ds
    C3_shorts[na, nd] = C3_short
    C3_longs[na, nd] = C3_long
    tofs[na, nd] = tof
    

# Prints the combination numver its calculating
print( '\nDeparture days: %i.'     % ds    )
print( 'Arrival days: %i.'         % as_   )
print( 'Total Combinations: %i.'   % total )




normed_departures = (et_departures - et_departures[0])/(3600.0 * 24.0)
normed_arrivals   = (et_arrivals  - et_arrivals[0])/(3600.0 * 24.0)


# Generate departure and arrival date grids
dep_mesh, arr_mesh = np.meshgrid( normed_departures, normed_arrivals )

# Create levels arrays
c3_levels = np.arange( 10, 50, 5)

# c3_levels_long = np.arange( 1900, 3900, 100)

vinf_levels = np.arange( 0, 15, 1)

tof_levels = np.arange( 100, 4000, 200)

# dv_levels  = np.arange( 3, 20, 0.5)




# Define linewdith
lw = 1.5

'''
Create The Plots
'''

fig, ax = plt.subplots(figsize = ( 15, 10 ))

c0 = ax.contour(dep_mesh, arr_mesh, C3_shorts, levels = c3_levels, colors = 'm', linewidths = lw)
plt.clabel( c0, fmt = '%i')

# c4 = ax.contour(dep_mesh, arr_mesh, tofs, levels = tof_levels, colors = 'red', linewidths = lw*0.6)
# plt.clabel( c4, fmt = '%i')



plt.plot( [ 0 ], [ 0 ], 'm' )
plt.plot( [ 0 ], [ 0 ], 'red' )


plt.legend(
  [
      'C3 ($\dfrac{km^2}{s^2}$)',
      'Time of Flight (days)'
  ],
  bbox_to_anchor = ( 1.005, 1.01 ),
  fontsize = 10
  )


# ax.set_ylim([0, 600])
ax.set_title('Porkchop Plot from Earth to Mars', fontsize = 20 )
ax.set_ylabel( 'Arrival (Days Past {})'.format(arrival0) , fontsize = 15 )
ax.set_xlabel( 'Departure (Days Past {})'.format(departure0), fontsize = 15 )
plt.show()
