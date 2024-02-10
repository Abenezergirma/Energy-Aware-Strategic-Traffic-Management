import matplotlib.pyplot as plt
import numpy as np
import json
from collections import defaultdict

plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42

capacity = 1
bar_width = 0.35
color1 = '#1f77b4'
color2 = '#ff7f0e'

fig, axs = plt.subplots(1, 3, figsize=(11, 4))
axs = axs.ravel()

bars = []

global_min_delay = float('inf')
global_max_delay = float('-inf')
global_min_energy = float('inf')
global_max_energy = float('-inf')
traffic = ['High','Medium','Low']

for idx, interval in enumerate([30, 60, 120]):
    with open(f'./ground_delay/delay_{interval}_{capacity}.json', 'r') as file:
        data = json.load(file)
    # Use zip to iterate over values in parallel based on the keys
    weights, ground_delays, energies = data['weight'], data['ground_delay'], data['energy']

    # Initialize defaultdicts for aggregation
    avg_ground_delay = defaultdict(list)
    avg_energy = defaultdict(list)
    # Aggregate data by weight
    for weight, ground_delay, energy in zip(weights, ground_delays, energies):
        avg_ground_delay[weight].append(ground_delay)
        avg_energy[weight].append(energy)

    avg_ground_delay = {w: sum(gd) / len(gd) for w, gd in avg_ground_delay.items()}
    avg_energy = {w: sum(e) / len(e) for w, e in avg_energy.items()}

    global_min_delay = min(0, min(avg_ground_delay.values())) #global_min_delay
    global_max_delay = max(850, max(avg_ground_delay.values())) #global_max_delay
    global_min_energy = min(0, min(avg_energy.values())) #global_min_energy
    global_max_energy = max(300, max(avg_energy.values())) #global_max_energy

    weights = sorted(avg_ground_delay.keys())
    avg_ground_delay_values = [avg_ground_delay[w] for w in weights]
    avg_energy_values = [avg_energy[w] for w in weights]

    index = np.arange(len(weights))
    ax = axs[idx]

    bar1 = ax.bar(index, avg_ground_delay_values, bar_width, label='Ground Delay', color=color1)
    ax.set_xlabel('Trade-off Weight')
    ax.set_ylabel(f'Avg Ground Delay (sec)', color=color1)
    ax.tick_params(axis='y', labelcolor=color1)
    ax.grid(True, linestyle='--', linewidth=0.5)
    ax.set_ylim(global_min_delay, global_max_delay)

    ax2 = ax.twinx()
    bar2 = ax2.bar(index + bar_width, avg_energy_values, bar_width, label='Energy Consumption', color=color2)
    ax2.set_ylabel('Avg Energy (MJ)', color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.set_ylim(global_min_energy, global_max_energy)

    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(weights)
    # Set a title for each subplot
    ax.set_title(f" {traffic[idx]} Traffic Density")

    if idx == 0:
        bars.extend([bar1, bar2])

# fig.suptitle(r'\textbf{Delay Vs Energy Tradeoff for Different Traffic Demands}', y=0.95)
fig.legend(handles=bars, labels=['Ground Delay', 'Energy Consumption'], loc='upper center', ncol=2, bbox_to_anchor=(0.5, 0.15))

plt.tight_layout(rect=[0, 0.1, 1, 0.95])
plt.savefig('combined_weight_results.pdf', bbox_inches='tight')
plt.show()
