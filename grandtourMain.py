# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 17:46:25 2025

@author: alexa
"""




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
step_size = 20

# --------------------

# Frame of reference from the planetary data
OBSERVER= pd.sun['name']
FRAME='ECLIPJ2000'
sun_mu = pd.sun['mu']



#cut off velocity (if the lambert solver fails to converge)
cutoff_v = 100
cutoff_c3 = cutoff_v ** 2

#Convert from days into seconds
step = step_size*3600*24



         

#Define function to get the emphemeris data - array holding positions and velocity vectors
def calc_ephemeris(target, ets, frame, observer):
	return np.array(spice.spkezr( target, ets, frame, 'NONE', observer )[ 0 ])

	
def norm(vec):
    return np.linalg.norm(vec)



def johann(dep_planet, arr_planet, departure0, departure1, arrival0, arrival1):
    
    et_departures = np.arange(spice.utc2et(departure0), spice.utc2et(departure1) + step, step)
    et_arrivals   = np.arange(spice.utc2et(arrival0), spice.utc2et(arrival1) + step, step)
    
    
    ds = len(et_departures)
    as_ = len(et_arrivals)
    C3_shorts     = np.zeros( (as_, ds) )
    C3_longs      = np.zeros( (as_, ds) )
    #tofs          = np.zeros( (as_, ds) )
    total = as_ * ds   
    
    # Retrieve the position and velocity vectors and store in a sixth-lengthed vector from function above
    ephem_departures = [calc_ephemeris(departure_planet, et, FRAME, OBSERVER) for et in et_departures]
    ephem_arrivals = [calc_ephemeris(arrival_planet, et, FRAME, OBSERVER) for et in et_arrivals]
    
    
    
    
    def compute_lambert_entry(na, nd):
        
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
        for na in tqdm(range(as_), desc="Calculating Transfer")
        for nd in range(ds))

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

minium value finder

'''


# min_index = np.unravel_index(np.argmin(johann_1[2]), johann_1[2].shape)
# min_value = johann_1[min_index]
# row, col = map(int, min_index)
# print(row, col)


# print(f"Miniumum energy occured at row {10*row} and column {10*col}")
# row_plus = 10*row + 370
# row_minus = 10*row - 70

# print(f"time for next trip is from {row_minus} and {row_plus}")

# date_stri = arrival_dates_0_i
# date_obji = datetime.strptime(date_stri, '%Y-%m-%d')

# # Add 70 days
# new_date_minus = date_obji + timedelta(days=row_minus)
# departure_dates_1_i = new_date_minus.strftime('%Y-%m-%d')

# new_date_plus = date_obji + timedelta(days=row_plus)
# departure_dates_1_f = new_date_plus.strftime('%Y-%m-%d')


# n_arrivaldate_lower = new_date_plus + timedelta(days=20)
# arrival_dates_1_i = n_arrivaldate_lower.strftime('%Y-%m-%d')

# n_arrivaldate_higher = new_date_plus + timedelta(days=1820)
# arrival_dates_1_f = n_arrivaldate_higher.strftime('%Y-%m-%d')

# print(departure_dates_1_i)
# print(departure_dates_1_f)
# print(arrival_dates_1_i)
# print(arrival_dates_1_f)


    
planet0 = 'Earth' #Case sensitive
planet1 = 'Mars'
# Dates for Launchwindow - wikipedia


dep_dates_0 = '2005-06-20'         # Intial departure date
dep_dates_1 = '2010-11-07'         # Final departure date
arr_dates_0 = '2005-12-01'         # Initial arrival date
arr_dates_1 = '2015-02-24'          # Final arrival date 
# Holds the spice name of the planets  
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
johann_1 = johann(departure_planet, arrival_planet, dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1)   


 
planet0 = 'Mars' #Case sensitive
planet1 = 'Jupiter'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_1_i = '2006-02-01'       
dep_dates_1_f = '2020-12-01'  
arr_dates_1_i = '2007-08-01'         
arr_dates_1_f = '2025-01-01'
johann_2 = johann(departure_planet, arrival_planet, dep_dates_1_i, dep_dates_1_f, arr_dates_1_i, arr_dates_1_f)   



planet0 = 'Jupiter' #Case sensitive
planet1 = 'Saturn'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_2_i =  "2009-01-01"
dep_dates_2_f = "2030-01-01"
arr_dates_2_i = "2011-01-01"
arr_dates_2_f = "2040-01-01"
johann_3 = johann(departure_planet, arrival_planet, dep_dates_2_i, dep_dates_2_f, arr_dates_2_i, arr_dates_2_f)   



date_axis = [(johann_1[0], johann_1[1]), (johann_2[0], johann_2[1]), (johann_3[0], johann_3[1])]
c3_shorts_arrays = [johann_1[2], johann_2[2], johann_3[2]]





'''
plot snake step 
'''


# Define linewdith
lw = 1.5
fig = plt.figure(figsize=(20,20))

n_plots = 3  # Adjust this as needed

w, h = 0.3, 0.3  # Size of each subplot
x, y = 0.1, 0.1  # Initial position

axes = []

c3_levels_0 = np.arange(10, 55, 10)
c3_levels_1 = np.arange( 10, 500, 50)
c3_levels_2 = np.arange( 10, 500, 100)

energy_levels = [c3_levels_0, c3_levels_1, c3_levels_2]


segment_names = ["Earth to Mars", "Mars to Jupiter", "Jupiter to Saturn"]
date_names = (
    (dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1),
    (dep_dates_1_i, dep_dates_1_f, arr_dates_1_i, arr_dates_1_f),
    (dep_dates_2_i, dep_dates_2_f, arr_dates_2_i, arr_dates_2_f)
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
        ax.set_xlabel(f"Departure (Days Past {date_names[i][0]})")
        ax.set_ylabel(f"Arrival (Days Past {date_names[i][2]})")
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.grid(True, linestyle='--', color='gray', linewidth=0.5)
        # if i == 0:
        #     plt.ylim(0, 400)
        #     ax.axhline(y=200, color='red', linestyle='--')
            
    # Odd plots: rotated
    else:
        ax.grid(True)
        # ref_ax = axes[i-1]
        # matching_yticks = ref_ax.get_yticks()
        Z_rot = np.rot90(c3_shorts_arrays[i])
        cont  = ax.contour(date_axis[i][1], date_axis[i][0][::-1],  Z_rot, levels=energy_levels[i],colors='m', linewidths = lw)
        plt.clabel( cont, fmt = '%i')
        
        ax.set_ylabel(f"Departure (Days Past {date_names[i][0]})")
        ax.set_xlabel(f"Arrival (Days Past {date_names[i][2]})")
        ax.xaxis.set_ticks_position('bottom')
        ax.xaxis.set_label_position('bottom')
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position('right')
        # if i == 1:
        #     plt.ylim(100, 300)
        #     ax.axhline(y=200, color='red', linestyle='--')
        #ax.set_yticks([])
        #ax.spines['left'].set_visible(False)
        
        

    ax.set_title(f"Segment {segment_names[i]}", fontsize=8)
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5)
    axes.append(ax)

    # Snake step layout
    if i % 2 == 0:
        x += w
    else:
        y += h


# for tick in matching_yticks:
#         # Convert from data to axes coordinates
#         ax.axhline(y=tick, linestyle='--', color='gray', linewidth=0.5, zorder=0)


# Optional: share gridlines
# for i in range(n_plots - 1):
#     if i % 2 == 0:
#         source_ax = axes[i]
#         target_ax = axes[i + 1]
#         for y_tick in np.linspace(0, 1, 6):  # Scaled y for visualization
#             target_ax.axhline(y=y_tick, linestyle='--', color='gray', linewidth=0.5, zorder=0)

#plt.suptitle("Snakestepping Porkchop Plot from Earth to Mars", fontsize=16)
plt.show()
