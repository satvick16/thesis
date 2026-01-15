import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Replace with your CSV file path
csv_file = "C:\\Users\\16474\\Desktop\\thesis\\touchstone\\pam3_driver_output.csv"

# Load CSV data
data = pd.read_csv(csv_file)

# Extract X and Y values
x = data['/vdrv_a X'].values
y = data['/vdrv_a Y'].values

# Compute the numerical derivatives
dy_dx = np.gradient(y, x)        # First derivative
d2y_dx2 = np.gradient(dy_dx, x) # Second derivative

# Create the plots
plt.figure(figsize=(10, 10))

# Original data
plt.subplot(3, 1, 1)
plt.plot(x, y, marker='.', linestyle='-', color='blue')
plt.title('Original Data')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)

# First derivative
plt.subplot(3, 1, 2)
plt.plot(x, dy_dx, marker='.', linestyle='-', color='red')
plt.title('First Derivative dY/dX')
plt.xlabel('X')
plt.ylabel('dY/dX')
plt.grid(True)

# Second derivative
plt.subplot(3, 1, 3)
plt.plot(x, d2y_dx2, marker='.', linestyle='-', color='green')
plt.title('Second Derivative d²Y/dX²')
plt.xlabel('X')
plt.ylabel('d²Y/dX²')
plt.grid(True)

plt.tight_layout()
plt.show()
