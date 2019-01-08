# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# How long does it take a neurotransimtter to diffuse across a synaptic cleft?
# While we can analytically solve for the first passage time with some
# knowledge of the physics, we can also simulate it as a simple one-dimensional
# random walk. We can model the synapse as a box with a reflective barrier
# (transmiting neuron) and an absorpitve barrier (receptor). In our model, the
# neurotransmitter can't travel pas teither barrier. If it is at teh reflective
# barrier, it can only travel forward one step size. If it is at the absorbing
# barrier, the walk is over!

# Before we can begin, we need to first define some parameters of our
# simulation.
synapse_length = 20  # in units of steps.
n_sim = 100000  # The number of simulations we'll do.
step_length = 1  # How far one molecule can travel in one time step

# Now we can actually perform the simulation. We'll first make an empty
# vector then jump into actually iterating for each simulation.
num_steps = np.zeros(n_sim)

# For each simulation...
for i in range(n_sim):

    # Set up a step counter.
    sim_steps = 0
    curr_pos = 0

    # So long as we haven't hit the receptor...
    while num_steps[i] == 0:

        # Make sure we haven't stepped beyond the transmitting neuron.
        if curr_pos < 0:
            curr_pos = 0
        elif curr_pos == synapse_length:
            num_steps[i] = sim_steps
            break  # Note that this stops the loope from being executed.
        else:
            # Now determine which way we will step.
            flip = np.random.rand()

            # This is assuming equal probability of steps.
            if flip > 0.5:
                curr_pos += step_length
            else:
                curr_pos -= step_length
            sim_steps += 1

# We can (and should) convert the number of steps to time, as that is what is
# physiologically meaningful. We can make an order-of-magnitude estimate that
# the diffusion constant ( D ) of a small organic molecule such as dopamine to
# be on the order of  ∼5cm2/s . We know that from Fick's law of diffusion, the
# time to diffuse a given distance is
#
#                               t ~ x^2 / 2D
#
# where  x is the distance of interest and  t is time. Note that the prefactor
# of 2 is for a one-dimensional random walk. In our simulation, the time to
# diffuse 1nm (a single step) is therefore
#
#                              t ~ 1 / 2D ~ 1μs.
#
# Since this is approximately 1µs, we can just use our num_steps vector as our
# time in units of µs. Let's plot the result of our situation and compute the
# confidence interval.

# Plot the distribution
plt.figure()
plt.hist(num_steps, bins=100)
plt.xlabel('time (µs)')
plt.ylabel('counts')
plt.show()

# Print the confidence interval and mean.
conf_int = np.percentile(num_steps, (2.5, 97.5))
mean_time = np.mean(num_steps)
print("""
The mean first-passage time is {0:.2f} µs.
The 2.5 percentile is {1:.2f} µs.
The 97.5 percentile is  {2:.2f} µs.
""".format(mean_time, conf_int[1], conf_int[0]))

# Notice that the diffusion across the synaptic cleft is much, much faster
# than an action potential!
