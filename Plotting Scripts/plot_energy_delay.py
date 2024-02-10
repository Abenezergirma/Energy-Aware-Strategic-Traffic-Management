import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import json
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42

# Load data (make sure to replace the file path with the correct one)
capacity = 1
interval = 60
with open(f'delay_{interval}_{capacity}.json', 'r') as file:
    data = json.load(file)

avg_ground_delay = defaultdict(list)
avg_energy = defaultdict(list)

for w, gd, e in zip(data['weight'], data['ground_delay'], data['energy']):
    # print(w)
    avg_ground_delay[w].append(gd)
    avg_energy[w].append(e)

avg_ground_delay = {w: sum(gd) / len(gd) for w, gd in avg_ground_delay.items()}
avg_energy = {w: sum(e) / len(e) for w, e in avg_energy.items()}

# Preparing data for plotting
weights = sorted(avg_ground_delay.keys())
avg_ground_delay_values = [avg_ground_delay[w] for w in weights]
avg_energy_values = [avg_energy[w] for w in weights]

# Plotting
bar_width = 0.35
index = np.arange(len(weights))
fig, ax1 = plt.subplots(figsize=(8, 4))  # Adjusted for double column
color1 = '#1f77b4'  # More sophisticated blue
color2 = '#ff7f0e'  # More sophisticated orange

# Plotting the first bar chart
ax1.bar(index, avg_ground_delay_values, bar_width, label='Ground Delay', color=color1)
ax1.set_xlabel('Weight')
ax1.set_ylabel('Average Ground Delay (seconds)', color=color1)
ax1.tick_params(axis='y', labelcolor=color1)
ax1.grid(True, linestyle='--', linewidth=0.5)

# Create a second y-axis for the second bar chart
ax2 = ax1.twinx()
ax2.bar(index + bar_width, avg_energy_values, bar_width, label='Energy Cost', color=color2)
ax2.set_ylabel('Average Energy (MJ)', color=color2)
ax2.tick_params(axis='y', labelcolor=color2)
ax2.grid(False)

# Common settings
plt.xticks(index + bar_width / 2, weights)
ax1.legend(loc='upper right', bbox_to_anchor=(1, 0.925))
ax2.legend(loc='upper right')

plt.tight_layout()  # Ensures good layout
plt.title('Delay Vs Energy Tradeoff')

# Save the figure
plt.savefig('weight_result.pdf', bbox_inches='tight')
plt.show()
