# Import our modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# Import our custom modules
import pboc_utils as pboc

# Set some parameter values.
r = 1  # Production rage in mRNA / min
gamma = 1 / 3  # in units of 1/min
time = 10  # total time of integration in min
dt = 0.05  # time step in min
num_steps = int(time / dt)
upper_bound = 20  # maximum mRNA copy number

# Make our storage array
prob = np.zeros((upper_bound + 1, num_steps))
prob[0, 0] = 1.0

# Spread the butter
for t in range(1, num_steps):
    for m in range(upper_bound):
        # Compute the probabilities
        prob[m, t] = prob[m, t-1] + gamma * (m + 1) * dt * prob[m+1, t-1] -\
                r * dt * prob[m, t-1] - gamma * m * dt * prob[m, t-1]

        # Consider when m > 0
        if m > 0:
            prob[m, t] = prob[m, t] + r * dt * prob[m - 1, t - 1]

# Change the time vector to real time.
time_vec = np.arange(0, time, dt)

pboc.bar3(prob, xlabel='time (min)', ylabel='number of mRNAs', zlabel='P(m,t)',
        x_vec=time_vec)
plt.show()
