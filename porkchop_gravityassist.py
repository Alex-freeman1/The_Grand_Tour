# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 12:14:36 2025

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
import planet_dates_ga as pg

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
    
    step_size = 100
    
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
planet1 = 'Jupiter'
# Holds the spice name of the planets  
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_0 = '1977-02-01'
dep_dates_1 = '1978-11-07'
arr_dates_0 = '1978-02-01'
arr_dates_1 = '1980-12-01'
johann_1 = johann(departure_planet, arrival_planet, dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1)   


'''
Jupiter - Saturn
''' 
planet0 = 'Jupiter' 
planet1 = 'Saturn'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
dep_dates_2_i = "1978-02-01"
dep_dates_2_f = "1980-12-01"
arr_dates_2_i = "1980-01-01"
arr_dates_2_f = "1983-01-01"
johann_3 = johann(departure_planet, arrival_planet, dep_dates_2_i, dep_dates_2_f, arr_dates_2_i, arr_dates_2_f)   


'''
Saturn - Uranus
''' 
planet0 = 'Saturn' 
planet1 = 'Uranus'
departure_planet, arrival_planet = loop_bodies(planet0, planet1)
arrival_planet = "URANUS BARYCENTER"
dep_dates_3_i =  "1980-01-01"
dep_dates_3_f = "1983-01-01"
arr_dates_3_i = "1986-01-01"
arr_dates_3_f = "1998-01-01"
johann_4 = johann(departure_planet, arrival_planet, dep_dates_3_i, dep_dates_3_f, arr_dates_3_i, arr_dates_3_f)   




# '''
# Uranus - Neptune
# ''' 
# departure_planet = "URANUS BARYCENTER"
# arrival_planet = "NEPTUNE BARYCENTER"
# dep_dates_4_i =  "1974-01-01"
# dep_dates_4_f = "2036-01-01"
# arr_dates_4_i = "1984-01-01"
# arr_dates_4_f = "2045-01-01"
# johann_5 = johann(departure_planet, arrival_planet, dep_dates_4_i, dep_dates_4_f, arr_dates_4_i, arr_dates_4_f)   




# Create arrays to hold the data from the johann function, all three indicies
date_axis = [(johann_1[0], johann_1[1]), (johann_3[0], johann_3[1]), (johann_4[0], johann_4[1])]
c3_shorts_arrays = [johann_1[2], johann_3[2], johann_4[2]]

'''
plot snake step 
'''


step_size__ = 50
start_date0 = spice.utc2et("1977-02-01")
start_date1 = spice.utc2et("1978-02-01")
start_date2 = spice.utc2et("1980-01-01")

ter_dates = pg.ter_dates_ga
jup_dates = pg.jup_dates_ga
sat_dates = pg.sat_dates_ga
#print(ter_dates[85],jup_dates[85], sat_dates[85])


points_0 = np.array(ter_dates) - start_date0
points_1 = np.array(jup_dates) - start_date1
points_2 = np.array(sat_dates) - start_date2


points_days_0 = points_0 / (86400)
points_days_1 = points_1 / (86400)
points_days_2 = points_2 / (86400)


# Define linewdith
lw = 0.5



# Let the number plots be a variable n_plots
n_plots = 2
fig = plt.figure(figsize=(100,100))   
# Size of each subplot
w, h = 0.4,0.4

# Initial position
x, y = 0.1, 0.15 

axes = []                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
font_size = 12
# Define levels for the contour plot 
c3_levels_earth = np.arange(0, 300, 20)
c3_levels_jupiter = np.arange(0, 300, 30)



energy_levels = [c3_levels_earth, c3_levels_jupiter, c3_levels_jupiter]

# Name the different segments
segment_names = ["Earth to Jupiter", "Jupiter to Saturn", "Saturn to Uranus", "Uranus to Neptune"]

# Create a tuple to hold all the dates for each porkchop plot
date_names = (
    (dep_dates_0, dep_dates_1, arr_dates_0, arr_dates_1),  
    (dep_dates_2_i, dep_dates_2_f, arr_dates_2_i, arr_dates_2_f),
    (dep_dates_3_i, dep_dates_3_f, arr_dates_3_i, arr_dates_3_f)
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
        x_end = ax.get_xlim()[1]
        y_end = ax.get_ylim()[1]
        
        t = [4,85]
        if i == 0:
            
            for n in t:
                color_x = 'red' if n == 4 else 'blue'
                ax.plot([points_days_0[n], x_end], [points_days_1[n], points_days_1[n]], color=color_x, linestyle='--')
                plt.plot(points_days_0[n], points_days_1[n], marker='^', color=color_x)
            
            # ax.axvline(x = points_days_0[0], color='red', linestyle='--')
            # ax.axhline(y = points_days_1[0], color='red', linestyle='--')
            
        if i == 2:
            
           for n in t:
                color_x = 'red' if n == 4 else 'blue'
                ax.axvline(x = points_days_2[n], color=color_x, linestyle='--')
            
            # for i in range(points_days_2.shape[0]):
            #     if points_days_2[i] < 0:
            #         pass
            #     else:
            #         ax.axvline(x = points_days_2[i], color='red', linestyle='--')
        
        ax.set_title(f"Segment {segment_names[i]}", fontsize=font_size)
        
    # Odd plots: rotated and flipped
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
        x_end = ax.get_xlim()[1]
        y_end = ax.get_ylim()[1]
        t = [4,85]
        if i == 1: 
            
            for n in t:
                color_x = 'red' if n == 4 else 'blue'
                ax.plot([0, points_days_2[n]], [points_days_1[n], points_days_1[n]], color=color_x, linestyle='--')
                ax.plot([points_days_2[n], points_days_2[n]], [points_days_1[n], y_end], color=color_x, linestyle='--')
                plt.plot(points_days_2[n], points_days_1[n], marker='^', color=color_x)
                
                
        ax.text(
            0.5, -0.1, f"Segment {segment_names[i]}",
            transform=ax.transAxes,
            ha='center', va='bottom', fontsize=font_size
        )

    # Snake step layout
    if i % 2 == 0:
        x += w
    else:
        y += h

plt.show()

