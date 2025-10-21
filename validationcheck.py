# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 17:46:25 2025

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


# Step size (in days)
#step_size = 50

# --------------------

# Frame of reference from the planetary data
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']



#cut off velocity (if the lambert solver fails to converge)
cutoff_v = 100
cutoff_c3 = cutoff_v ** 2

#Convert from days into seconds



#Define function to get the emphemeris data - array holding positions and velocity vectors
def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

	
# Returns the normalised vector
def norm(vec):
    return np.linalg.norm(vec)

# Function that uses parallel computing to loop through each planet journey and calculate the lambert transfers
# It returns the C_3 energy of each transfer in a 2 dimensional array with the axis as the departure and arrival dates
def johann(dep_planet, arr_planet, departure0, departure1, arrival0, arrival1):
      
    step_size = 2
    
    step = step_size*3600*24
    
    
    # Sets the depature and arrival arrays using the SPICE utc2et function
    et_departures = np.arange(spice.utc2et(departure0), spice.utc2et(departure1) + step, step)
    et_arrivals   = np.arange(spice.utc2et(arrival0), spice.utc2et(arrival1) + step, step)
    
    # Gets length of these arrays
    ds = len(et_departures)
    as_ = len(et_arrivals)
    
    # Begins 2 dimensional arrays to hold the data
    C3_shorts     = np.zeros( (as_, ds) )
    C3_longs      = np.zeros( (as_, ds) )
    #tofs          = np.zeros( (as_, ds) )
    
    # Calculates how many entries exist
    total = as_ * ds   
    
    # Retrieve the position and velocity vectors and store in a sixth-lengthed vector from function above
    ephem_departures = [calc_ephemeris(departure_planet, et, FRAME, OBSERVER) for et in et_departures]
    ephem_arrivals = [calc_ephemeris(arrival_planet, et, FRAME, OBSERVER) for et in et_arrivals]
    
    # Function to pass into the parallel computing job function
    def compute_lambert_entry(states_depart, states_arrive, tof, sun_mu, cutoff_c3):
        from numpy.linalg import norm
        
        
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
        delayed(compute_lambert_entry)(
            ephem_departures[nd],
            ephem_arrivals[na],
            et_arrivals[na] - et_departures[nd],
            sun_mu,
            cutoff_c3
        )
        for na in tqdm(range(as_), desc="Calculating Transfer")
        for nd in range(ds)
)

    for idx, (C3_short, C3_long, tof) in enumerate(results):
        na = idx // ds
        nd = idx % ds
        C3_shorts[na, nd] = C3_short
        C3_longs[na, nd] = C3_long
        #tofs[na, nd] = tof

    
    # Prints the combination numver its calculating
    print( '\nDeparture days: %i.'     % ds    )
    print( 'Arrival days: %i.'         % as_   )
    print( 'Total Combinations: %i.'   % total )
    
    # Convert tof from sec to days
    #tofs /= ( 3600.0 * 24.0 )
    
    # Total delta-v
    # dv_shorts = v_inf_shorts + np.sqrt( C3_shorts )
    # dv_longs  = v_inf_longs  + np.sqrt( C3_longs  )
    
    
    normed_departures = (et_departures - et_departures[0])/(3600.0 * 24.0)
    normed_arrivals   = (et_arrivals  - et_arrivals[0])/(3600.0 * 24.0)

    return normed_departures, normed_arrivals, C3_shorts 
 


# Loops through the podies to check against the inputted planet0 and planet1
pd_bodies = pd.bodies
def loop_bodies(planeti, planetf):
    pd_req0 = pd_req1 = None
    for body in pd_bodies:
        if body['name'] == planeti:
            pd_req0 = body['spice_name']
        elif body['name'] == planetf:
            pd_req1 = body['spice_name']
        
    return pd_req0, pd_req1
        
   
 


'''
Earth - Mars
'''    
planet0 = 'Earth' #Case sensitive
planet1 = 'Mars'
# Holds the spice name of the planets  
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_0 = '2005-06-20'         # Intial departure date
dep_dates_1 = '2005-11-07'         # Final departure date
arr_dates_0 = '2005-12-01'         # Initial arrival date
arr_dates_1 = '2007-02-24'          # Final arrival date 
johann_1 = johann(departure_planet, arrival_planet, dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1)   



# Create arrays to hold the data from the johann function, all three indicies

johanns = [johann_1]

date_axis = []
c3_shorts_arrays = []

for j in johanns:
    date_axis.append((j[0], j[1]))
    c3_shorts_arrays.append(j[2])
    


'''
plot snake step 
'''

# Define linewdith
lw = 0.3
fig = plt.figure(figsize=(100,100))   


# Let the number plots be a variable n_plots
n_plots = len(johanns)

# Size of each subplot
w, h = 0.4,0.8

# Initial position
x, y = 0.1, 0.1

axes = []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
font_size = 15
# Define levels for the contour plot 

c3_levels_0 = np.arange(0, 30, 2)



energy_levels = [c3_levels_0]

# Name the different segments
segment_names = ["Earth to Mars"]
# Create a tuple to hold all the dates for each porkchop plot
date_names = (
    (dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1),  

)
for i in range(n_plots):
    
    ax = fig.add_axes([x, y, w, h])
    plt.legend(
      [
          'C3 ($\dfrac{km^2}{s^2}$)',
         
      ],
      bbox_to_anchor = ( 1.005, 1.01 ),
      fontsize = 10
      )
    # Even plots: original orientation
    if i % 2 == 0:
        
        cont  = ax.contour(date_axis[i][0], date_axis[i][1], c3_shorts_arrays[i], levels=energy_levels[i], colors='m', linewidths = lw)
        plt.clabel( cont, fmt = '%i')
        ax.set_xlabel(f"Departure (Days Past {date_names[i][0]})", fontsize=font_size)
        ax.set_ylabel(f"Arrival (Days Past {date_names[i][2]})", fontsize=font_size)
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.grid(True, linestyle='--', color='gray', linewidth=0.5)
        
        ax.set_title(f"Segment {segment_names[i]}", fontsize=font_size)

    # Odd plots: rotated and flipped
    else:
        ax.grid(True)
        # ref_ax = axes[i-1]
        # matching_yticks = ref_ax.get_yticks()
        Z_rot = np.rot90(c3_shorts_arrays[i])
        cont  = ax.contour(date_axis[i][1], date_axis[i][0][::-1],  Z_rot, levels=energy_levels[i],colors='m', linewidths = lw)
        plt.clabel( cont, fmt = '%i')
        
        ax.set_ylabel(f"Departure (Days Past {date_names[i][0]})", fontsize=font_size)
        ax.set_xlabel(f"Arrival (Days Past {date_names[i][2]})", fontsize=font_size)
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_label_position('bottom')
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')
    
        
    
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5)
    axes.append(ax) 
    plt.tick_params(axis='both', which='major', labelsize=10)

    # Snake step layout
    if i % 2 == 0:
        x += w
    else:
        y += h

plt.show()
