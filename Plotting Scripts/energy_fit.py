import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42


# Define the datasets
data = pd.DataFrame({
    'Speed': [46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
    'Total_Optimal_Cost': [352.653, 337.834, 323.905, 310.865, 298.576, 286.998, 276.083, 265.776, 256.039, 246.826, 238.086, 229.831, 221.984]
})

data2 = pd.DataFrame({
    'Speed': [46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
    'Total_Optimal_Cost': [352.696, 337.876, 323.969, 310.903, 298.613, 287.034, 276.116, 265.81, 256.071, 246.856, 238.116, 229.859, 222.011]
})

data3 = pd.DataFrame({
    'Speed': [46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58],
    'Total_Optimal_Cost': [249.967, 239.463, 229.608, 220.346, 211.635, 203.43, 195.694, 188.388, 181.485, 174.954, 168.773, 162.908, 157.346]
})

# Function to fit model and return predictions
def fit_and_predict(data):
    model = LinearRegression()
    x = data['Speed'].values.reshape(-1, 1)
    y = data['Total_Optimal_Cost'].values.reshape(-1, 1)
    model.fit(x, y)
    predictions = model.predict(x)
    return x, y, predictions

# Fit models and get predictions
x1, y1, predictions1 = fit_and_predict(data)
x2, y2, predictions2 = fit_and_predict(data2)
x3, y3, predictions3 = fit_and_predict(data3)

# Plotting
plt.figure(figsize=(4, 3)) 
plt.scatter(x1, y1, color='blue', label='Raw Data 1',marker='o')
plt.plot(x1, predictions1, color='red', label='Fitted Curve 1')

plt.scatter(x2, y2, color='green', label='Raw Data 2',marker='^')
plt.plot(x2, predictions2, color='orange', label='Fitted Curve 2')

plt.scatter(x3, y3, color='purple', label='Raw Data 3',marker='s')
plt.plot(x3, predictions3, color='brown', label='Fitted Curve 3')

plt.xlabel('Speed (m/s)', fontsize=10)
plt.ylabel('Energy consumption (MJ)', fontsize=10)
# plt.title(r'\textbf{Energy Vs Speed Comparison}', fontsize=12)
plt.legend(loc='upper right',fontsize=8)
plt.grid(True, linestyle='--', alpha=0.7)
# Save the figure
plt.savefig('energy_speed.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()
