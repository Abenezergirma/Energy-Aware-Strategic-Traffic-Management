import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
plt.rcParams['pdf.fonttype'] = 42

# Define the range of x, avoiding zero to prevent division by zero
x = np.linspace(5, 70)
x2 = np.linspace(40, 60)
x3 = np.linspace(30, 40)
y = 1 / x
x_point = 50
x_point2 = 35

a1 = 1/x_point
a2 = 1/(x_point*x_point)
y_taylor = a1 - a2 * (x2 - x_point)
y_point = a1

a3 = 1/x_point2
a4 = 1/(x_point2*x_point2)
y_taylor2 = a3 - a4 * (x3 - x_point2)
y_point2 = a3

# Plotting the curve
plt.figure(figsize=(4, 3))
plt.plot(x, y, label='$y = 1/x$')
plt.scatter(x_point, y_point, color='red')
# plt.scatter(x_point2, y_point2, color='red')
plt.plot(x2, y_taylor, label='Taylor series expension')
# plt.plot(x3, y_taylor2, label='Taylor series expension')
plt.title('Taylor Series Expansion of $y = 1/x$')
plt.xlabel('x')
plt.ylabel('y')
# plt.grid(True)
plt.grid(True, linestyle='--', linewidth=0.5)
plt.legend()
# Save the figure
plt.savefig('taylor_series.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()
