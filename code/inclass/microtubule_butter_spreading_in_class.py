# ...
import numpy as np
import matplotlib.pyplot as plt
import seaborn
import pboc_utils as pboc

# Define the parameters.
r = 5  # addition rate in monomers per time
gamma = 6  # Degradation rate in monomers per time
dt = 1 / 50 # Time step for spreading.
tot_length = 100  # Maximum polymer length
tot_time = 1000   # Time to integrate in time steps

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
#pboc.bar3(prob, xlabel='time (steps)', ylabel='polymer length in monomers',\
#        zlabel = 'P(l, t)', bin_step=50)
#plt.show()

# Compute the mean length at the end of the integration.
lengths = np.arange(0, tot_length + 1, 1)
length_prob = lengths * prob[:, -1]
mean_length = np.sum(length_prob)
predicted_length = (r / gamma) / (1 - (r / gamma))
print('Our mean length is ' + str(mean_length) + '. Our predicted length is ' + str(predicted_length))

# Show the distribution at the end.
plt.figure()
plt.bar(lengths, prob[:, -1])

# Our predicted distribution.
prediction = (1 - (r / gamma)) * (r / gamma)**lengths
plt.plot(lengths, prediction, 'r-')
plt.xlabel('length in monomers')
plt.ylabel('P(l)')
plt.show()


# Now consider a length dependent gamma.

# Set up the two-dimensional matrix of probabilities.
prob = np.zeros((tot_length + 1, tot_time))

# Set our initial condition.
prob[30, 0] = 1.0

# Set a constant for our length dependent gamma.
c = 1 / 10
r = 1

# Start the integration.
for t in range(1, tot_time):
    for ell in range(tot_length):
        # Evaluate the master equation.
        prob[ell, t] = prob[ell, t-1] - r * dt * prob[ell, t-1] +\
            c * (ell + 1) * dt * prob[ell+1, t-1]

        if ell > 0:
            # prob[ell, t] += x is same as prob[ell,  t] = prob[ell, t] + x
            prob[ell, t] = prob[ell, t] + r * dt * prob[ell-1, t-1] - c * ell * dt * prob[ell, t-1]

# Show the distribution.
pboc.bar3(prob, xlabel='time (steps)', ylabel='length in monomers', zlabel='$P(\ell)$',
        bin_step=50)
plt.show()
