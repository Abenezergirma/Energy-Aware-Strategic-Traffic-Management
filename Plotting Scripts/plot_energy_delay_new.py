import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import json

plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42

capacity = 1
bar_width = 0.35
color1 = '#1f77b4'  # More sophisticated blue
color2 = '#ff7f0e'  # More sophisticated orange

# Define the font size
axis_label_fontsize = 14
tick_label_fontsize = 12

# Create a 2x2 subplot layout
fig, axs = plt.subplots(2, 2, figsize=(12, 10))  # Adjust the figure size as needed
axs = axs.ravel()  # Flatten the 2D array of axes for easy indexing

bars = []  # To store bar objects for the legend

# Loop through each interval and plot
for idx, interval in enumerate([30, 60, 120, 240]):
    with open(f'./ground_delay/delay_{interval}_{capacity}.json', 'r') as file:
        data = json.load(file)

    avg_ground_delay = defaultdict(list)
    avg_energy = defaultdict(list)

    for w, gd, e in zip(data['weight'], data['ground_delay'], data['energy']):
        avg_ground_delay[w].append(gd)
        avg_energy[w].append(e)

    avg_ground_delay = {w: sum(gd) / len(gd) for w, gd in avg_ground_delay.items()}
    avg_energy = {w: sum(e) / len(e) for w, e in avg_energy.items()}

    weights = sorted(avg_ground_delay.keys())
    avg_ground_delay_values = [avg_ground_delay[w] for w in weights]
    avg_energy_values = [avg_energy[w] for w in weights]

    index = np.arange(len(weights))
    ax = axs[idx]  # Select the subplot

    # Plotting the first bar chart
    bar1 = ax.bar(index, avg_ground_delay_values, bar_width, label='Ground Delay', color=color1)
    ax.set_xlabel('Weight', fontsize=axis_label_fontsize)
    ax.set_ylabel(f'Avg Ground Delay (sec)', color=color1, fontsize=axis_label_fontsize)
    ax.tick_params(axis='both', labelsize=tick_label_fontsize)
    ax.grid(True, linestyle='--', linewidth=0.5)

    # Create a second y-axis for the second bar chart
    ax2 = ax.twinx()
    bar2 = ax2.bar(index + bar_width, avg_energy_values, bar_width, label='Energy Cost', color=color2)
    ax2.set_ylabel('Avg Energy (MJ)', color=color2, fontsize=axis_label_fontsize)
    ax2.tick_params(axis='both', labelsize=tick_label_fontsize)

    ax.set_xticks(index + bar_width / 2, weights)

    if idx == 0:  # Only add to the legend list once
        bars.extend([bar1, bar2])

fig.suptitle(r'\textbf{Delay Vs Energy Tradeoff for Different Intervals}', fontsize=16, fontweight='bold')


# Place a single shared legend outside of the subplots
fig.legend(handles=bars, labels=['Ground Delay', 'Energy Cost'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.08), fontsize=14)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout to fit
plt.subplots_adjust(top=0.95)  # Adjust the top value as needed
plt.savefig('combined_weight_results.pdf', bbox_inches='tight')
plt.show()

