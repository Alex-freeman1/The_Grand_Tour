


'''
Lambert's Problem Tools
'''

# Third-party Libraries
import numpy as np

'''
This is helpful Python script to run necessary functions to solve Lambert's Problem.
The specific details for this sciprt were found in the report done by Saad and Nouh found:
    
'A.N. Saad and M.I. Nouh, (2003). Lambert Universal Variable Algorithm, 
Arabian Journal for Science and Engineering, 88-94.'

This derives all the variables that this script uses
# 
'''

# This is a function which is derived from Newton's Interation Method,
# It takes in some x value as an initial guess and changes it slightly until some 
# Condition is met; in this case it's when the x value is 
# Within the tolerance
def newton_interation(time, dydx, init_guess, tol, i=1000):
    x = init_guess

    for i in range(i):
        x_new = x - time(x) / dydx(x)
        if np.abs(x_new - x) < tol:
            return x_new
        x = x_new
    


# Complete Universal Lambert Solver as laid out in the research paper
# It takes in the position vectors with respect to some frame, in this case it will be the sun,
# The time of flight, as well as the Gravitational potential factor, mu, again in this case it 
# will be the sun
def lambert_solver(R1, R2, dt, mu, tol=1e-6, maxiter=10000, trajectory='pro'):

    
# Define the second and third stumphff functions that are used in celestial orbital mechanics 
# As per the Wikipedia page
    def stumphff2(x):
        if x > 0:
            return (1 - np.cos(np.sqrt(x))) / x
        elif x < 0:
            return (np.cosh(np.sqrt(-x)) - 1) / (-x)
        else:
            return 1 / 2
    def stumphff3(x):
        if x > 0:
            return (np.sqrt(x) - np.sin(np.sqrt(x))) / (np.sqrt(x)) ** 3
        elif x < 0:
            return (np.sinh(np.sqrt(-x)) - np.sqrt(-x)) / (np.sqrt(-x)) ** 3
        else:
            return 1 / 6

    # Define the coeffecient B as used in the paper
    def coeff_B(x):
        return r1 + r2 + A * (x * stumphff3(x) - 1) / np.sqrt(stumphff2(x))

    # Define the function to return the change in time of flight
    def delta_T(x):
        return (coeff_B(x) / stumphff2(x)) ** 1.5 * stumphff3(x) + A *  np.sqrt(coeff_B(x)) - np.sqrt(mu) * dt
    
    # Function to pass into the intetation method 
    def change_tol(x):
        if x == 0:
            return np.sqrt(2) / 40 * coeff_B(0) ** 1.5 + A / 8 * (np.sqrt(coeff_B(0)) + A * np.sqrt(1 / 2 / coeff_B(0)))
        else:
            return (coeff_B(x) / stumphff2(x)) ** 1.5 * (1 / 2 / x * (stumphff2(x) - 3 * stumphff3(x) / 2 / stumphff2(x)) + 3 * stumphff3(x) ** 2 / 4 / stumphff2(x)) + A / 8 * (3 * stumphff3(x) / stumphff2(x) * np.sqrt(coeff_B(x) ) + A * np.sqrt(stumphff2(x) / coeff_B(x) ))
        


    # Magnitudes of R1 and R2
    r1 = np.linalg.norm(R1)
    r2 = np.linalg.norm(R2)

    
    # Find the cross product of these two R-3 vectors
    cross12 = np.cross(R1, R2)
    
    # Calculate the angle between the vectors and normalise it
    dtheta = np.arccos(np.dot(R1, R2) / (r1 * r2))
    
    
    # Depending on whether we are taking the prograde or retrograde path,
    # We need to manipulate the angle if it is negative or positive respectively, 
    # So that it is in some range between 0 and 2pi
    if trajectory == 'pro':  # For prograde trajectory
        if cross12[2] < 0:
            dtheta = 2 * np.pi - dtheta
    elif trajectory == 'retro':  # For retrograde trajectory
        if cross12[2] >= 0:
            dtheta = 2 * np.pi - dtheta
    else:
        pass

    # Calulate A as per in the research paper
    A = np.sin(dtheta) * np.sqrt(r1 * r2 / (1 - np.cos(dtheta)))
    
    
    # Initialise some variable x to compute the Newton Interation
    x = 0.1 
    while delta_T(x) < 0:
        x += 0.1
        

    # Complete the Newton Interation function
    x = newton_interation(delta_T, change_tol, x, tol, maxiter)


    # Check if x is some reasonable number
    if np.isnan(x) or np.isinf(x):
        solved = False
    else:
        solved = True 
    
    
    
    if solved:
        
        f = 1 - coeff_B(x) / r1
        #fdot = (np.sqrt(mu) / (r1 * r2)) * np.sqrt(coeff_B(x) / stumphff2(x)) * (x * stumphff3(x) - 1)
        g = A * np.sqrt(coeff_B(x) / mu)
        gdot = 1 - coeff_B(x) / r2
    
        # Compute the velocities V1 & V2
        V1 = 1 / g * (R2 - f * R1)
        V2 = 1 / g * (gdot * R2 - R1)
            
        # Return the velocity vectors for leaving Earth and arriving 
        # At Mars respectively
        return V1, V2
    
