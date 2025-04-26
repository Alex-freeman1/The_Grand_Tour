# -*- coding: utf-8 -*-
"""
Created on Wed Apr 23 17:58:53 2025

@author: alexa
"""
import matplotlib.pyplot as plt
import numpy as np

# Example data
fig = plt.figure(figsize=(10, 10))
size = 5
Z = np.random.rand(size, size)
departure_dates = np.linspace(0, 10, size)
arrival_dates = np.linspace(0, 20, size)



n_plots = 5
w, h = 0.3, 0.3  # Size of each subplot
x, y = 0.1, 0.1  # Starting position

axes = []

for i in range(n_plots):
    ax = fig.add_axes([x, y, w, h])
    
    # Alternate between original and rotated view
    if i % 2 == 0:
        ax.imshow(Z, extent=[departure_dates[0], departure_dates[-1],
                             arrival_dates[0], arrival_dates[-1]],
                  aspect='auto', origin='lower')
        ax.set_xlabel("Departure")
        ax.set_ylabel("Arrival")
    else:
        Z_rot = np.rot90(Z)
        ax.imshow(Z_rot, extent=[arrival_dates[0], arrival_dates[-1],
                                 departure_dates[-1], departure_dates[0]],
                  aspect='auto', origin='lower')
        ax.set_xlabel("Arrival")
        ax.set_ylabel("Departure")
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_yticks([])  # Hide right y-ticks
        ax.spines['left'].set_visible(False)

    ax.set_title(f"Plot {i+1}", fontsize=8)
    ax.grid(True, linestyle='--', color='gray', linewidth=0.5)

    axes.append(ax)

    # Step snake layout
    if i % 2 == 0:
        x += w
    else:
        y += h

# Extend horizontal grid lines from left to right plots
y_ticks = np.linspace(0, 1, 6)
for i in range(n_plots - 1):
    if i % 2 == 0:
        source_ax = axes[i]
        target_ax = axes[i + 1]
        for y_tick in y_ticks:
            target_ax.axhline(y=y_tick, linestyle='--', color='gray', linewidth=0.5, zorder=0)

plt.show()
