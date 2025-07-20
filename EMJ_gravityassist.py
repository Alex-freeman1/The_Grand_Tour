
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import spiceypy as spice
import lambertsolve_master as lt
import planetary_data as pd
from scipy.integrate import solve_ivp

# Load SPICE kernels
spice.furnsh("de432s\de432s.bsp")
spice.furnsh("data\latest_leapseconds (1).tls")

def seconds(days):
    return days * 24 * 3600

def days(seconds_val):
    return seconds_val / (24 * 3600)

def calc_ephemeris(target, ets, frame, observer):
    return np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])

def norm(vec):
    return np.linalg.norm(vec)

def calc_vinfinity(tof, target_planet, departure_planet, et0, vinf_target):
    """
    Function to compute difference between target vinfinity magnitude and
    current vinfinity magnitude for a given time of flight (tof).
    """
    r1_target = calc_ephemeris(target_planet, et0 + tof, FRAME, OBSERVER)[:3]
    state_departure = calc_ephemeris(departure_planet, et0, FRAME, OBSERVER)
    
    v0_sc_depart, v1_sc_arrive = lt.lambert_solver(
        state_departure[:3],    
        r1_target,   
        tof,
        sun_mu,
        trajectory='pro'
    )
    
    vinf = np.linalg.norm(v0_sc_depart - state_departure[3:])
    return vinf_target - vinf

def vinfinity_match(planet0, planet1, v_sc_incoming, et0, tof0, diff_step=1e-3, tol=1e-4):
    """
    Match V-infinity magnitudes for hyperbolic flyby
    """
    state0_planet0 = calc_ephemeris(planet0, et0, FRAME, OBSERVER)
    vinf_incoming = v_sc_incoming - state0_planet0[3:]
    vinf_magnitude = np.linalg.norm(vinf_incoming)
    
   
    def root_func(tof):
        return calc_vinfinity(tof, planet1, planet0, et0, vinf_magnitude)
    
    tof, steps = lt.newton_root_single_fd(root_func, tof0, tol=tol, diff_step=diff_step)
    
    r1_planet1 = calc_ephemeris(planet1, et0 + tof, FRAME, OBSERVER)[:3]
    
    v_sc_depart, v_sc_arrive = lt.lambert_solver(
        state0_planet0[:3],    
        r1_planet1,   
        tof,
        sun_mu,
        trajectory='pro'
    )
    
    return tof, v_sc_depart, v_sc_arrive



# Constants
OBSERVER = pd.sun['name']
FRAME = 'ECLIPJ2000'
sun_mu = pd.sun['mu']

# Get planetary data
pd_bodies = pd.bodies
for i in range(10):
    if pd_bodies[i]['name'] == 'Earth':
        pd_earth = pd_bodies[i]
    elif pd_bodies[i]['name'] == 'Jupiter':
        pd_mars = pd_bodies[i]
    elif pd_bodies[i]['name'] == 'Saturn':
        pd_jupiter = pd_bodies[i]

Earth = pd_earth['spice_name']
Jupiter = pd_mars['spice_name']
Saturn = pd_jupiter['spice_name']

# # EMJ trajectory times 
# time0 = "2029-10-01"  # Earth departure
# time1 = "2030-08-15"  # Mars flyby
# time2 = "2034-07-07"  # Jupiter arrival

#EMJ trajectory times 

time0 = "1977-08-15"  # Earth departure
time1 = "1979-07-22"  # Jupiter flyby
time2 = "1981-08-07"  # Saturn arrival

et_0 = spice.utc2et(time0)
et_1 = spice.utc2et(time1)
et_2 = spice.utc2et(time2)

print("EMJ Trajectory Analysis")


# LEG 1: Earth to Mars
print("LEG 1: Earth to Mars")
state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
state1 = calc_ephemeris(Jupiter, et_1, FRAME, OBSERVER)
tof_1 = et_1 - et_0

print(f"  Time of flight: {days(tof_1):.1f} days")

# Lambert solution for Earth-Jupiter leg
v_sc_depart_1, v_sc_arrive_1 = lt.lambert_solver(
    state0[:3],    
    state1[:3],   
    tof_1,
    sun_mu,
    trajectory='pro'
)

# Earth departure analysis
v_inf_earth_departure = v_sc_depart_1 - state0[3:]


#print(f"  Earth departure V-infinity: {norm(v_inf_earth_departure):.3f} km/s")


# Jupiter arrival analysis
v_inf_mars_arrival = v_sc_arrive_1 - state1[3:]
print(f"Jupiter arrival V-infinity: {norm(v_inf_mars_arrival):.3f} km/s")
print()


#  ----------------------------------------------------------------------------

# LEG 2: Jupiter to Saturn 
print("LEG 2: Jupiter to Saturn")
tof_guess_2 = et_2 - et_1  # Initial guess

try:
    tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
        Jupiter, Saturn,
        v_sc_arrive_1,  # spacecraft velocity at Mars arrival
        et_1, tof_guess_2
    )
    
    print(f"Converged! Time of flight: {days(tof_2):.1f} days")
    
    # Verify the Jupiter flyby
    state_Jupiter_flyby = calc_ephemeris(Jupiter, et_1, FRAME, OBSERVER)
    vinf_in = v_sc_arrive_1 - state_Jupiter_flyby[3:]
    vinf_out = v_sc_depart_2 - state_Jupiter_flyby[3:]
    
    print(v_sc_depart_2)
    #print(state_mars_flyby[3:])
    print("\n")
    
    # print(f"Jupiter flyby incoming V-inf: {vinf_in} km/s")
    # print(f"Jupiter flyby outgoing V-inf: {vinf_out} km/s")
    
    # Calculate deflection angle
    deflection_angle = np.arccos(np.dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))
    print(f"Mars deflection angle: {np.degrees(deflection_angle):.1f} degrees")
    
    # Saturn arrival analysis
    state_Saturn_arrival = calc_ephemeris(Saturn, et_1 + tof_2, FRAME, OBSERVER)
    v_inf_Saturn_arrival = v_sc_arrive_2 - state_Saturn_arrival[3:]
    
    
    

except Exception as e:
    print(f"  Error in Mars-Jupiter leg: {e}")
    print("  Trying alternative approaches...")
    
    # Try different time of flight guesses
    for alt_days in [500, 800, 1000, 1200]:
        print(f"  Trying {alt_days} day transfer...")
        try:
            tof_2_alt = seconds(alt_days)
            tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
                Jupiter, Saturn,
                v_sc_arrive_1,
                et_1, tof_2_alt
            )
            print(f"  SUCCESS with {alt_days} day guess: actual TOF = {days(tof_2):.1f} days")
            break
        except Exception as e2:
            print(f"    Failed: {e2}")
    else:
        print("  All alternative guesses failed. Consider:")
        print("  1. Different departure/arrival dates")
        print("  2. Different trajectory type (retrograde)")
        print("  3. Multiple Mars flybys")



v_inf = norm(vinf_in)
print(vinf_in)
print(vinf_out)
print(v_sc_arrive_1)
print(v_sc_depart_2)

r_jup = 69911
mu_jup = 126.687*10**6

def altitude(delta_a, v_infx):
    e = 1 / np.sin(delta_a / 2)
    r_p = (e - 1) * mu_jup / v_inf**2
    h = r_p - r_jup
    return h

h = altitude(np.radians(60.9), v_inf)
print(f"Flyby altitude above Mars surface: {h} km")


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Mars at origin
ax.quiver(0, 0, 0, *vinf_in, color='red', label='v_inf_in')
ax.quiver(0, 0, 0, *vinf_out, color='blue', label='v_inf_out')

# Axis limits (adjust as needed)
ax.set_xlim([-20, 20])
ax.set_ylim([-20, 20])
ax.set_zlim([-10, 10])

ax.set_xlabel('Vx (km/s)')
ax.set_ylabel('Vy (km/s)')
ax.set_zlabel('Vz (km/s)')
ax.set_title('Hyperbolic Excess Velocity Vectors at Mars')
ax.legend()
plt.show()



# def rotation_matrix_z(theta):
#     return np.array([
#         [np.cos(theta), -np.sin(theta), 0],
#         [np.sin(theta),  np.cos(theta), 0],
#         [0,              0,             1]
#     ])




# print("\n")
# print(vinf_in)
# print(vinf_out)
# R_z = rotation_matrix_z(deflection_angle)

# print(vinf_in)
# print(deflection_angle)
# v_inf_plus = R_z @ vinf_in
# print(v_inf_plus)






# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# import spiceypy as spice
# import lambertsolve_master as lt
# import planetary_data as pd
# from scipy.integrate import solve_ivp

# # Load SPICE kernels
# spice.furnsh("de432s\de432s.bsp")
# spice.furnsh("data\latest_leapseconds (1).tls")

# def seconds(days):
#     return days * 24 * 3600

# def days(seconds_val):
#     return seconds_val / (24 * 3600)

# def calc_ephemeris(target, ets, frame, observer):
#     return np.array(spice.spkezr(target, ets, frame, 'NONE', observer)[0])

# def norm(vec):
#     return np.linalg.norm(vec)

# def calc_vinfinity(tof, target_planet, departure_planet, et0, vinf_target):
#     """
#     Function to compute difference between target vinfinity magnitude and
#     current vinfinity magnitude for a given time of flight (tof).
#     """
#     r1_target = calc_ephemeris(target_planet, et0 + tof, FRAME, OBSERVER)[:3]
#     state_departure = calc_ephemeris(departure_planet, et0, FRAME, OBSERVER)
    
#     v0_sc_depart, v1_sc_arrive = lt.lambert_solver(
#         state_departure[:3],    
#         r1_target,   
#         tof,
#         sun_mu,
#         trajectory='pro'
#     )
    
#     vinf = np.linalg.norm(v0_sc_depart - state_departure[3:])
#     return vinf_target - vinf

# def vinfinity_match(planet0, planet1, v_sc_incoming, et0, tof0, diff_step=1e-3, tol=1e-4):
#     """
#     Match V-infinity magnitudes for hyperbolic flyby
#     """
#     state0_planet0 = calc_ephemeris(planet0, et0, FRAME, OBSERVER)
#     vinf_incoming = v_sc_incoming - state0_planet0[3:]
#     vinf_magnitude = np.linalg.norm(vinf_incoming)
    
   
#     def root_func(tof):
#         return calc_vinfinity(tof, planet1, planet0, et0, vinf_magnitude)
    
#     tof, steps = lt.newton_root_single_fd(root_func, tof0, tol=tol, diff_step=diff_step)
    
#     r1_planet1 = calc_ephemeris(planet1, et0 + tof, FRAME, OBSERVER)[:3]
    
#     v_sc_depart, v_sc_arrive = lt.lambert_solver(
#         state0_planet0[:3],    
#         r1_planet1,   
#         tof,
#         sun_mu,
#         trajectory='pro'
#     )
    
#     return tof, v_sc_depart, v_sc_arrive



# # Constants
# OBSERVER = pd.sun['name']
# FRAME = 'ECLIPJ2000'
# sun_mu = pd.sun['mu']

# # Get planetary data
# pd_bodies = pd.bodies
# for i in range(10):
#     if pd_bodies[i]['name'] == 'Earth':
#         pd_earth = pd_bodies[i]
#     elif pd_bodies[i]['name'] == 'Mars':
#         pd_mars = pd_bodies[i]
#     elif pd_bodies[i]['name'] == 'Jupiter':
#         pd_jupiter = pd_bodies[i]

# Earth = pd_earth['spice_name']
# Mars = pd_mars['spice_name']
# Jupiter = pd_jupiter['spice_name']

# # # EMJ trajectory times 
# # time0 = "2029-10-01"  # Earth departure
# # time1 = "2030-08-15"  # Mars flyby
# # time2 = "2034-07-07"  # Jupiter arrival

# #EMJ trajectory times 

# time0 = "2007-12-18"  # Earth departure
# time1 = "2008-06-22"  # Mars flyby
# time2 = "2013-01-07"  # Jupiter arrival

# et_0 = spice.utc2et(time0)
# et_1 = spice.utc2et(time1)
# et_2 = spice.utc2et(time2)

# print("EMJ Trajectory Analysis")


# # LEG 1: Earth to Mars
# print("LEG 1: Earth to Mars")
# state0 = calc_ephemeris(Earth, et_0, FRAME, OBSERVER)
# state1 = calc_ephemeris(Mars, et_1, FRAME, OBSERVER)
# tof_1 = et_1 - et_0

# print(f"  Time of flight: {days(tof_1):.1f} days")

# # Lambert solution for Earth-Mars leg
# v_sc_depart_1, v_sc_arrive_1 = lt.lambert_solver(
#     state0[:3],    
#     state1[:3],   
#     tof_1,
#     sun_mu,
#     trajectory='pro'
# )

# # Earth departure analysis
# v_inf_earth_departure = v_sc_depart_1 - state0[3:]


# #print(f"  Earth departure V-infinity: {norm(v_inf_earth_departure):.3f} km/s")


# # Mars arrival analysis
# v_inf_mars_arrival = v_sc_arrive_1 - state1[3:]
# print(f"  Mars arrival V-infinity: {norm(v_inf_mars_arrival):.3f} km/s")
# print()


# #  ----------------------------------------------------------------------------

# # LEG 2: Mars to Jupiter 
# print("LEG 2: Mars to Jupiter")
# tof_guess_2 = et_2 - et_1  # Initial guess

# try:
#     tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
#         Mars, Jupiter,
#         v_sc_arrive_1,  # spacecraft velocity at Mars arrival
#         et_1, tof_guess_2
#     )
    
#     print(f"  Converged! Time of flight: {days(tof_2):.1f} days")
    
#     # Verify the Mars flyby
#     state_mars_flyby = calc_ephemeris(Mars, et_1, FRAME, OBSERVER)
#     vinf_in = v_sc_arrive_1 - state_mars_flyby[3:]
#     vinf_out = v_sc_depart_2 - state_mars_flyby[3:]
    
#     print(v_sc_depart_2)
#     #print(state_mars_flyby[3:])
#     print("\n")
    
#     # print(f"Mars flyby incoming V-inf: {vinf_in} km/s")
#     # print(f"Mars flyby outgoing V-inf: {vinf_out} km/s")
    
#     # Calculate deflection angle
#     deflection_angle = np.arccos(np.dot(vinf_in, vinf_out) / (norm(vinf_in) * norm(vinf_out)))
#     print(f"Mars deflection angle: {np.degrees(deflection_angle):.1f} degrees")
    
#     # Jupiter arrival analysis
#     state_jupiter_arrival = calc_ephemeris(Jupiter, et_1 + tof_2, FRAME, OBSERVER)
#     v_inf_jupiter_arrival = v_sc_arrive_2 - state_jupiter_arrival[3:]
    
    
    

# except Exception as e:
#     print(f"  Error in Mars-Jupiter leg: {e}")
#     print("  Trying alternative approaches...")
    
#     # Try different time of flight guesses
#     for alt_days in [500, 800, 1000, 1200]:
#         print(f"  Trying {alt_days} day transfer...")
#         try:
#             tof_2_alt = seconds(alt_days)
#             tof_2, v_sc_depart_2, v_sc_arrive_2 = vinfinity_match(
#                 Mars, Jupiter,
#                 v_sc_arrive_1,
#                 et_1, tof_2_alt
#             )
#             print(f"  SUCCESS with {alt_days} day guess: actual TOF = {days(tof_2):.1f} days")
#             break
#         except Exception as e2:
#             print(f"    Failed: {e2}")
#     else:
#         print("  All alternative guesses failed. Consider:")
#         print("  1. Different departure/arrival dates")
#         print("  2. Different trajectory type (retrograde)")
#         print("  3. Multiple Mars flybys")



# v_inf = norm(vinf_in)
# print(v_inf)
# v_infplus = norm(vinf_out)
# print(v_infplus)

# print(v_sc_arrive_1)
# print(v_sc_depart_2)

# r_mars = 3389.5 
# mu_mars = 42828.3

# def altitude(delta_a, v_infx):
#     e = 1 / np.sin(delta_a / 2)
#     r_p = (e - 1) * mu_mars / v_inf**2
#     h = r_p - r_mars
#     return h

# h = altitude(np.radians(60.9), v_inf)
# print(f"Flyby altitude above Mars surface: {h:.2f} km")


# # fig = plt.figure()
# # ax = fig.add_subplot(111, projection='3d')

# # # Mars at origin
# # ax.quiver(0, 0, 0, *vinf_in, color='red', label='v_inf_in')
# # ax.quiver(0, 0, 0, *vinf_out, color='blue', label='v_inf_out')

# # # Axis limits (adjust as needed)
# # ax.set_xlim([-20, 20])
# # ax.set_ylim([-20, 20])
# # ax.set_zlim([-10, 10])

# # ax.set_xlabel('Vx (km/s)')
# # ax.set_ylabel('Vy (km/s)')
# # ax.set_zlabel('Vz (km/s)')
# # ax.set_title('Hyperbolic Excess Velocity Vectors at Mars')
# # ax.legend()
# # plt.show()



# # def rotation_matrix_z(theta):
# #     return np.array([
# #         [np.cos(theta), -np.sin(theta), 0],
# #         [np.sin(theta),  np.cos(theta), 0],
# #         [0,              0,             1]
# #     ])




# # print("\n")
# # print(vinf_in)
# # print(vinf_out)
# # R_z = rotation_matrix_z(deflection_angle)

# # print(vinf_in)
# # print(deflection_angle)
# # v_inf_plus = R_z @ vinf_in
# # print(v_inf_plus)










     





     
