import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42
# Define the path to your .dat file
file_path = 'u_no_windalt_500speed_55.dat'#'x_no_windalt_500speed_55.dat'


data = []  # List to hold all arrays

try:
    with open(file_path, 'r') as file:
        for line in file:
            # Assuming the values are tab-separated
            values = line.strip().split('\t')
            
            # Convert the values from strings to floats
            float_values = [float(value) for value in values]
            print(float_values)
            data.append(float_values)

except FileNotFoundError:
    print(f"File not found: {file_path}")
except ValueError:
    print("Encountered a value that couldn't be converted to float.")
except Exception as e:
    print(f"An error occurred: {e}")

# Check if we have at least three rows
if len(data) >= 3:
    # Plot the second and third rows
    plt.figure(figsize=(10, 6))
    plt.plot(data[1], label='Second Row')
    plt.plot(data[2], label='Third Row')

    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.title('Plot of 2nd and 3rd Rows from File')
    plt.legend()
    plt.show()
else:
    print("Not enough rows in the data file.")
