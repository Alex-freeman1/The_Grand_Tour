import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(10, 10))

n_plots = 5
w, h = 0.3, 0.3  # Size of each subplot

x, y = 0.1, 0.1
axes = []

# Set fixed y-limits and y-ticks
y_min, y_max = 0, 1
y_ticks = np.linspace(y_min, y_max, 6)

for i in range(n_plots):
    ax = fig.add_axes([x, y, w, h])
    ax.plot([0, 1], [0, 1])
    
    ax.set_ylim(y_min, y_max)
    ax.set_yticks(y_ticks)
    ax.grid(True, which='both', linestyle='--', color='gray', linewidth=0.5)
    ax.set_title(f"{i+1}", fontsize=8)
    
    axes.append(ax)

    if i % 2 == 0:
        x += w
    else:
        y += h
        ax.xaxis.set_ticks_position('top')
        ax.xaxis.set_label_position('top')
        ax.set_yticks([])  # Hide y-ticks for right-hand axes
        ax.spines['left'].set_visible(False)

# Extend horizontal grid lines from odd to even plots
for i in range(n_plots - 1):
    if i % 2 == 0:
        source_ax = axes[i]
        target_ax = axes[i + 1]

        for y_tick in y_ticks:
            target_ax.axhline(y=y_tick, linestyle='--', color='gray', linewidth=0.5, zorder=0)

plt.show()
