# Import our modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn
# Import our custom module
import pboc_utils as pboc

# Define the parameter values.
k = 5  # rate in  1 / s
dt = 1 / 50   # time step in seconds

# Set up the array for the probabilities.
num_boxes = 100
time_points = 100
prob = np.zeros((num_boxes, time_points))

# Set our initial condition
prob[49, 0] = 1.0

# Evaluate the master equation
for t in range(1, time_points):
    for x in range(1, num_boxes - 1):
        prob[x, t] = prob[x, t - 1] + k * dt * prob[x - 1, t - 1] -\
         2 * k * dt * prob[x, t - 1] + k * dt * prob[x + 1, t - 1]

pboc.bar3(prob, xlabel='time (s)', ylabel='box number', zlabel='probability',
          bin_step=5)

# Set up the array for the probabilities.
num_boxes = 10
time_points = 100
prob = np.zeros((num_boxes, time_points))

# Set our initial condition
prob[:, 0] = 1.0 / 4
prob[2:8, 0] = 0

# Evaluate the master equation
for t in range(1, time_points):
    # Set the boundary conditions
    prob[0, t] = prob[0, t - 1] - k * dt * prob[0, t - 1] +\
                k * dt * prob[1, t - 1]

    # Set condition for the last box
    prob[-1, t] = prob[-1, t - 1] + k * dt * prob[-2, t - 1] - k * dt *prob[-1, t-1]

    for x in range(1, num_boxes - 1):
        prob[x, t] = prob[x, t - 1] + k * dt * prob[x - 1, t - 1] -\
         2 * k * dt * prob[x, t - 1] + k * dt * prob[x + 1, t - 1]

# Plot it
pboc.bar3(prob, xlabel='time', ylabel='box number', zlabel='probability')

plt.show()
