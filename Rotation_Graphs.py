import matplotlib.pyplot as plt
import numpy as np

# Example data
Z = np.random.rand(2, 2)
departure_dates = np.linspace(0, 10, 10)
arrival_dates = np.linspace(0, 20, 10)

# Plot 1: Original
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.imshow(Z, extent=[departure_dates[0], departure_dates[-1],
                      arrival_dates[0], arrival_dates[-1]],
           aspect='auto', origin='lower')
plt.xlabel("Departure Date")
plt.ylabel("Arrival Date")
plt.title("Original")

# Plot 2: Rotated 90 degrees counterclockwise
plt.subplot(1, 2, 2)
Z_rotated = np.rot90(Z)  # Rotate the matrix 90°
plt.imshow(Z_rotated, extent=[arrival_dates[0], arrival_dates[-1],
                               departure_dates[-1], departure_dates[0]],
           aspect='auto', origin='lower')
plt.xlabel("Arrival Date")
plt.ylabel("Departure Date")
plt.title("Rotated 90°")

plt.tight_layout()
plt.show()
