# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pboc_utils as pboc

# In this script, we will look at the distribution of polymer ell_range under
# some different types of regulation. We will first think of some cytoskeletal
# polymer growing in an unlimited pool of available monomers (tubulin dimers)
# and the case in which the shortening of the polymer is dependent on the total
# length.

# To begin, we'll define a set of parameters. We'll assume that the production
# rate and the depoymerization rate are independent of length. Our master
# equation becomes
#
#  P(l, t + dt) = P(l, t) + r * dt * P(l-1, t) + y * dt * P(l+1, t) - r * dt * P(l, t + dt) - gamma * dt * P(l, t + dt)
#
# where l is the length of the polymer, t is the time, dt is the timestep, r
# is the growth rate of the polymer, and y is the depolymerization rate.
r = 5  # addition rate in monomers per time
gamma = 6  # Degradation rate in monomers per time
dt = 1 / 50 # Time step for spreading.
tot_length = 100  # Maximum polymer length
tot_time = 1500   # Time to integrate in time steps

# Now we are ready to perform integration as we typically do.

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

        # Now we can work out the boundary condition when ell is 0.
        if ell > 0:
            prob[ell, t] = prob[ell, t] + r * dt * prob[ell-1, t-1] - \
                            gamma * dt * prob[ell, t-1]

# Now we'll show the plot in three dimensions.
pboc.bar3(prob, xlabel='time (steps)', ylabel='polymer length in monomers',
          zlabel='P(l, t)', bin_step=50)
plt.show()

# We can see that we are approaching steady state and the distribution appears
# to be exponential. What is the mean?  From our master equation, we can reach
# an equation that describes the mean polymer length at steady state,
#
#   <ell> = (r/y)  / (1 - (r/y))
#
# We can compute the mean length of our polymers and compare it to the
# prediction.
ell_range = np.arange(0, tot_length + 1, 1)
ell_probabilities = ell_range * prob[:, -1]
mean_length = np.sum(ell_probabilities)
predicted_length = (r / gamma) / (1 - (r/gamma))
print('Our predicted length is ' + str(predicted_length) + '. Our computed mean length is ' + str(mean_length) + '.')


# We are actually pretty close, but not exactly correct. This could be due to
# numerical precision or maybe we aren't quite at setady state yet. However,
# this is pretty close. Let's now see how well the entire distribution agrees
# with our numerical integration. Remember, the entire distribution can be
# described by
#
#  P(l) = (r / y)**l * (1 - (r / y)).
#
predicted_distribution = (r / gamma)**ell_range * (1 - (r / gamma))
plt.figure()
plt.bar(ell_range, prob[:,-1])
plt.plot(ell_range, predicted_distribution, 'r-')
plt.xlabel('length in momonmers')
plt.ylabel('$P(l)$')
plt.show()

# That agrees quite well!


# Now, let's look at the case in which gamma is dependent on the length of the
# polymer. In essence, this means that we now can redfine gamma as some
# constant times the length of the polymer.
gamma_ell = 1/2
r = 1

# We'll reset our probability matrix and redo the integration.
prob = np.zeros((tot_length + 1, tot_time))

# We'll set our initial condition this time with long filaments.
prob[70, 0] = 1.0

# Now we just integrate.
for t in range(1, tot_time):
    for ell in range(tot_length):
        # Evaluate the master equation.
        prob[ell, t] = prob[ell, t-1] - r * dt * prob[ell, t-1] +\
            gamma_ell * (ell + 1) * dt * prob[ell+1, t-1]

        # Now we can work out the boundary condition when ell is 0.
        if ell > 0:
            prob[ell, t] = prob[ell, t] + r * dt * prob[ell-1, t-1] - \
                            gamma_ell * ell * dt * prob[ell, t-1]

# Now we'll show the plot in three dimensions.
pboc.bar3(prob, xlabel='time (steps)', ylabel='polymer length in monomers',
          zlabel='P(l, t)', bin_step=50)
plt.show()

# We can see that it now decreases in length to teh steady state value. As
# we did before, we can compare our computed mean length to the predicted
# value. For this model, the predicted mean length is
#
# <ell> = (r / gamma_ell)
#
length_prob = ell_range * prob[:, -1]
mean_length = np.sum(length_prob)
predicted_mean =  r / gamma_ell
print('Our predicted length is ' + str(predicted_mean) + '. Our computed length is ' + str(mean_length) + '.')


# Again, that agrees quite well. We can solve for the steady state distribution
# of polymer ell_range as
#
# P(l) = (r / y')^l * (1 / l!) * e^(-(r/y'))
#
predicted_dist = []
for ell in ell_range:
    prediction = (r / gamma_ell)**ell * (1 / np.math.factorial(ell)) *\
                  np.exp(-r/gamma_ell)
    predicted_dist.append(prediction)

# Let's plot it.
plt.figure()
plt.bar(ell_range, prob[:, -1])
plt.plot(ell_range, predicted_dist, 'r-')
plt.xlabel('polymer length in monomers')
plt.ylabel('$P(\ell)$')
plt.show()

# Again, it looks like we've nailed it.
