# ...
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set_context('talk')
import pboc_utils as pboc

# Define the parameters.
r = 19 # addition rate in monomers per time
gamma = 20  # Degradation rate in monomers per time
dt = 1 / 100  # Time step for spreading.
tot_length = 100  # Maximum polymer length
tot_time = 800   # Time to integrate in time steps

# Set up the two-dimensional matrix of probabilities.
prob = np.zeros((tot_length + 1, tot_time))

# Set our initial condition.
prob[0, 0] = 1.0

# Start the integration.
for t in range(1, tot_time):
    for ell in range(tot_length):
        # Evaluate the master equation.
        prob[ell, t] = prob[ell, t-1] - r * dt * prob[ell, t-1] +\
            gamma * dt * prob[ell+1, t-1]

        if ell > 0:
            # prob[ell, t] += x is same as prob[ell,  t] = prob[ell, t] + x
            prob[ell, t] = prob[ell, t] + r * dt * prob[ell-1, t-1] - gamma * dt * prob[ell, t-1]

# Generate the plot.
pboc.bar3(prob, xlabel='time (steps)', ylabel='polymer length in monomers',\
        zlabel = '$P(\ell, t)$', bin_step=10)
plt.show()
plt.tight_layout()
plt.savefig('../figures/microtubule_master_equation.png', bbox_inches='tight')
