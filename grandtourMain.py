# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 17:46:25 2025

@author: alexa
"""




#Import Python files
import matplotlib.pyplot as plt
import spiceypy as spice
import numpy as np

#Import Lambert tools
import lambertSolver as lt
import planetary_data as pd

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
departure0 = '2028-09-15'         
departure1 = '2029-06-01'  
arrival0 = '2029-05-01'         
arrival1 = '2030-04-01'

# Step size (in days)
step_size = 10

# --------------------

# Dates for Launchwindow - wikipedia
departure0 = '2005-06-20'         # Intial departure date
departure1 = '2005-11-07'         # Final departure date
arrival0 = '2005-12-01'         # Initial arrival date
arrival1 = '2008-02-24'          # Final arrival date

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

           
# For each combination of departure and arrivals
for na in y:
    for nd in x:
    
        
        # Retrieve the position and velocity vectors and store in a sixth-lengthed vector from function above
        states_depart = calc_ephemeris(departure_planet, et_departures[nd], FRAME, OBSERVER)
        states_arrive = calc_ephemeris(arrival_planet, et_arrivals[na], FRAME, OBSERVER)
        
        
        
        # The time of flight for each orbital transfer
        tof = (et_arrivals[na] - et_departures[nd])

        # Attempt to solve Lambert's problem for velocities

        # Short way (prograde)
        try:
            v_sc_depart_short, v_sc_arrive_short = lt.lambert_solver(
                # states_arrive is a 6 element vector, the first three represent the position vector: x,y,z
                # the final three represent the velocity vector: x,y,z.
                # therefore :3 is taking the first three position vectors passing it as R1
                states_depart[ :3 ],
                states_arrive[ :3 ], #R2
                tof,
                sun_mu,
                trajectory='pro'
            )
        except:
            v_sc_depart_short = np.array( [1000, 1000, 1000] )
            v_sc_arrive_short = np.array( [1000, 1000, 1000] )
            
            
        
        
        # Long way (retrograde)
        try:
            v_sc_depart_long, v_sc_arrive_long = lt.lambert_solver(
                states_depart[ :3 ],
                states_arrive[ :3 ],
                tof,
                sun_mu,
                trajectory='retro'
            )
        except:
           
            v_sc_depart_long = np.array( [1000, 1000, 1000] )
            v_sc_arrive_long = np.array( [1000, 1000, 1000] )
            
            
        
        
        # Calculate C3 values at departure
        
        # states_arrive is a 6 element vector, the first three represent the position vector: x,y,z
        # the final three represent the velocity vector: x,y,z.
        C3_short = norm(v_sc_depart_short - states_depart[ 3: ]) ** 2
        C3_long  = norm( v_sc_depart_long  - states_depart[ 3: ] ) ** 2

        # Check for unreasonable values (C3)
        if C3_short > cutoff_c3: 
            C3_short = cutoff_c3
        if C3_long  > cutoff_c3: 
            C3_long = cutoff_c3

        # Calculate v_infinity values at arrival
        v_inf_short = norm( v_sc_arrive_short - states_arrive[ 3: ] ) 
        v_inf_long  = norm( v_sc_arrive_long  - states_arrive[ 3: ] )

        # Check for unreasonable values (v_infinity)
        if v_inf_short > cutoff_v: 
            v_inf_short = cutoff_v
        if v_inf_long  > cutoff_v: 
            v_inf_long  = cutoff_v

        # Append values to corresponding arrays
        C3_shorts    [ na, nd ] = C3_short
        C3_longs     [ na, nd ] = C3_long
        # v_inf_shorts [ na, nd ] = v_inf_short
        # v_inf_longs  [ na, nd ] = v_inf_long
        tofs         [ na, nd ] = tof
       

    print( f'{na + 1} / {as_}.' )


# Prints the combination numver its calculating
print( '\nDeparture days: %i.'     % ds    )
print( 'Arrival days: %i.'         % as_   )
print( 'Total Combinations: %i.'   % total )

# Convert tof from sec to days
tofs /= ( 3600.0 * 24.0 )

# Total delta-v
# dv_shorts = v_inf_shorts + np.sqrt( C3_shorts )
# dv_longs  = v_inf_longs  + np.sqrt( C3_longs  )


normed_departures = (et_departures - et_departures[0])/(3600.0 * 24.0)
normed_arrivals   = (et_arrivals  - et_arrivals[0])/(3600.0 * 24.0)


# Generate departure and arrival date grids
dep_mesh, arr_mesh = np.meshgrid( normed_departures, normed_arrivals )

# Create levels arrays
c3_levels = np.arange( 10, 50, 2)

# c3_levels_long = np.arange( 1900, 3900, 100)

vinf_levels = np.arange( 0, 15, 1)

tof_levels = np.arange( 100, 600, 20)

# dv_levels  = np.arange( 3, 20, 0.5)







'''
break
'''

# Define linewdith
lw = 1.5
fig = plt.figure(figsize=(12, 12))

n_plots = 5  # Adjust this as needed
w, h = 0.3, 0.3  # Size of each subplot
x, y = 0.1, 0.1  # Initial position

axes = []
c3_levels = np.arange( 10, 50, 2)


for i in range(n_plots):
    ax = fig.add_axes([x, y, w, h])
    
    # Even plots: original orientation
    if i % 2 == 0:
        # Can use cmap='viridis' or equivalent to highlight edge cases if resolution faulters
        cont  = ax.contour(normed_departures, normed_arrivals, C3_shorts, levels=c3_levels, colors='m', linewidths = lw)
        ax.set_xlabel("Departure (Days)")
        ax.set_ylabel("Arrival (Days)")
        ax.grid(True, linestyle='--', color='gray', linewidth=0.5)
    # Odd plots: rotated
    else:
        ax.grid(True)
        # ref_ax = axes[i-1]
        # matching_yticks = ref_ax.get_yticks()
        Z_rot = np.rot90(C3_shorts)
        cont  = ax.contour(normed_arrivals, normed_departures[::-1],  Z_rot, levels=c3_levels,colors='m', linewidths = lw)
        ax.set_xlabel("Arrival (Days)")
        ax.set_ylabel("Departure (Days)")
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_yticks([])
        ax.spines['left'].set_visible(False)

    ax.set_title(f"Segment {i+1}", fontsize=8)
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

plt.suptitle("Snakestepping Porkchop Plot from Earth to Mars", fontsize=16)
plt.show()
