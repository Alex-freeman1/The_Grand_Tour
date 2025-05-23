


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
def newton_interation(time, dydx, init_guess, tol, i=100):
    x = init_guess

    for _ in range(i):
        fx = time(x)
        dfx = dydx(x)

        # Guard: avoid division by zero or nan/inf
        if dfx == 0 or np.isnan(dfx) or np.isinf(dfx):
            #print(f"Newton failed: derivative invalid at x = {x}")
            return np.nan

        x_new = x - fx / dfx

        if np.abs(x_new - x) < tol:
            return x_new

        x = x_new

    #print("Newton iteration did not converge.")
    return np.nan
    


# Complete Universal Lambert Solver as laid out in the research paper
# It takes in the position vectors with respect to some frame, in this case it will be the sun,
# The time of flight, as well as the Gravitational potential factor, mu, again in this case it 
# will be the sun
def lambert_solver(R1, R2, dt, mu, tol=1e-6, maxiter=1000, trajectory='pro'):

    
# Define the second and third stumphff functions that are used in celestial orbital mechanics 
# As per the Wikipedia page
    def stumphff2(x):
        eps = 1e-6
        threshold = 80  # Avoid overflow in cosh/sinh
    
        if x > eps:
            sqrt_x = np.sqrt(x)
            return (1 - np.cos(sqrt_x)) / x
        elif x < -eps:
            sqrt_neg_x = np.sqrt(-x)
            if sqrt_neg_x > threshold:
                return np.inf  # or return np.nan to indicate failure
            return (np.cosh(sqrt_neg_x) - 1) / (-x)
        else:
            # Use Taylor series around 0
            return 1 / 2 - x / 24 + x**2 / 720
    
    def stumphff3(x):
        eps = 1e-8
        threshold = 100  # Avoid overflow in sinh
    
        if x > eps:
            sqrt_x = np.sqrt(x)
            return (sqrt_x - np.sin(sqrt_x)) / (sqrt_x ** 3)
        elif x < -eps:
            sqrt_neg_x = np.sqrt(-x)
            if sqrt_neg_x > threshold:
                return np.inf  # or np.nan
            return (np.sinh(sqrt_neg_x) - sqrt_neg_x) / (sqrt_neg_x ** 3)
        else:
            return 1 / 6 - x / 120 + x**2 / 5040

    # Define the coeffecient B as used in the paper
    def coeff_B(x):
        s2 = stumphff2(x)
        if s2 <= 0 or np.isnan(s2) or np.isinf(s2):
            return np.nan  # or return some large number or failure code
        try:
            return r1 + r2 + A * (x * stumphff3(x) - 1) / np.sqrt(s2)
        except:
            return np.nan

    # Define the function to return the change in time of flight
    def delta_T(x):
        try:
            B = coeff_B(x)
            S = stumphff2(x)
            T = stumphff3(x)
    
            if S <= 0 or B <= 0 or not np.isfinite(B) or not np.isfinite(S):
                return np.inf  # Unphysical or unstable
    
            term1 = (B / S) ** 1.5 * T
            term2 = A * np.sqrt(B)
    
            F = term1 + term2
            if not np.isfinite(F):
                return np.inf
    
            return F - np.sqrt(mu) * dt
        except Exception:
            return np.inf
    
    
    # def delta_T(x):
    #     return (coeff_B(x) / stumphff2(x)) ** 1.5 * stumphff3(x) + A *  np.sqrt(coeff_B(x)) - np.sqrt(mu) * dt
    
    # Function to pass into the intetation method 
    def change_tol(x):
        try:
            B = coeff_B(x)
            S2 = stumphff2(x)
            S3 = stumphff3(x)
    
            # Prevent invalid or unsafe operations
            if x == 0 or S2 <= 0 or B <= 0 or not np.isfinite(B) or not np.isfinite(S2):
                return np.inf  # Or return 0 if you prefer a graceful fallback
    
            ratio = B / S2
            sqrtB = np.sqrt(B)
            sqrt_ratio = ratio ** 1.5
    
            # Construct derivative safely
            term1 = sqrt_ratio * (1 / (2 * x) * (S2 - 1.5 * S3 / S2) + 0.75 * S3**2 / S2)
            term2 = A / 8 * (3 * S3 / S2 * sqrtB + A * np.sqrt(S2 / B))
    
            df = term1 + term2
    
            return df if np.isfinite(df) else np.inf
        except Exception:
            return np.inf

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
    attempts = 0
    while delta_T(x) < 0:
        x += 0.1
        attempts += 1
        if attempts > 1000:
            print("No feasible Lambert solution: delta_T(x) < 0 for all x tried.")
            return None
        

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
    