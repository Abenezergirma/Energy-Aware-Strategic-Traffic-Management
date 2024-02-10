import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import json
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42

# Load data (make sure to replace the file path with the correct one)
capacity = 1
interval = 120
with open(f'new_profile_{interval}_{capacity}.json', 'r') as file:
    json_data = json.load(file)
    

# Preparing the data for plotting
data_by_weight = {}  # Dictionary to hold data by weight

for entry in json_data:
    weight = entry['weight']
    for profile in entry['profile']:
        speed = profile['speed']
        delay = profile['ground_delay']
        if weight not in data_by_weight:
            data_by_weight[weight] = {'speeds': [], 'delays': []}
        data_by_weight[weight]['speeds'].append(speed)
        data_by_weight[weight]['delays'].append(delay)

# Plotting
# plt.figure(figsize=(10, 6))

# for weight, data in data_by_weight.items():
    # if weight != 10:
    #     continue
    # plt.scatter(data['speeds'], data['delays'], label=f'Weight: {weight}')
    # print('speed ',data['speeds'])
    # print('delays ',data['delays'])



# Data
speed = [56, 53, 56, 50, 56, 56, 50, 56, 50, 47, 50, 50, 47, 50, 50, 50, 50, 50, 50, 50, 47, 50, 50, 47]

# Create a histogram of speeds
# plt.figure(figsize=(10, 6))
# plt.hist(speed, bins=range(min(speed), max(speed) + 2, 1), alpha=0.7, color='blue', edgecolor='black')

# # Customize the plot
# plt.xlabel('Speed')
# plt.ylabel('Number of Aircraft')
# plt.title('Distribution of Aircraft Speeds')
# plt.xticks(range(min(speed), max(speed) + 1))  # Ensures every speed is marked
# plt.grid(axis='y', alpha=0.75)

# # Show plot
# plt.show()

import matplotlib.pyplot as plt
import numpy as np

# Provided data
speed = [56, 53, 56, 50, 56, 56, 50, 56, 50, 47, 50, 50, 47, 50, 50, 50, 50, 50, 50, 50, 47, 50, 50, 47]
delays = [0, 0, 0, 10, 72, 71, 242, 80, 333, 292, 535, 356, 466, 443, 803, 692, 760, 695, 1036, 1182, 1357, 1083, 1374, 1312]

# Normalize the delay values to [0, 2*pi] for plotting
max_delay = max(delays)
theta = [2 * np.pi * delay / max_delay for delay in delays]

# Use speed as radial distance
r = speed

# Create polar plot
# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# ax.scatter(theta, r)

# # Add labels and title
# ax.set_title('Speed vs. Ground Delay')
# ax.set_theta_zero_location('N')  # Set the direction of 0 degrees to the top
# ax.set_theta_direction(-1)  # Set the direction of angles to clockwise
# ax.set_rlabel_position(90)  # Set the position of the radial labels

# plt.show()

# X locations for the groups
ind = np.arange(len(speed))  # the x locations for the groups
width = 0.35  # the width of the bars

# plt.figure(figsize=(4, 3)) 
fig, ax1 = plt.subplots(figsize=(4, 3))

# Bar plot for speed
rects1 = ax1.bar(ind - width/2, speed, width, label='Speed', color='SkyBlue')

# Create another y-axis for the delays
ax2 = ax1.twinx()
rects2 = ax2.bar(ind + width/2, delays, width, label='Delay', color='IndianRed')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax1.set_xlabel('Index of Aircraft')
ax1.set_ylabel('Speed', color='SkyBlue')
ax2.set_ylabel('Delay', color='IndianRed')
# ax1.set_title('Speed and Ground Delay Comparison')
ax1.set_xticks(ind)
ax1.set_xticklabels([str(i+1) for i in range(len(speed))], rotation=70)
ax1.tick_params(axis='x', labelsize=7) # Reduces the font size of x-axis tick labels
ax1.grid(True, linestyle='--', linewidth=0.5)
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')

ax1.tick_params(axis='y', labelcolor='SkyBlue')
ax2.tick_params(axis='y', labelcolor='IndianRed')

fig.tight_layout()  # To make sure everything fits without overlapping
plt.savefig('speed_delay.pdf', bbox_inches='tight')

plt.show()






