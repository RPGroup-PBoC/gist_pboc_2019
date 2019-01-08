# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# This script will perform a sotcastic random walk with a varying degree
# of bias in teh probability of taking a step. In our previous simulations,
# we have asserted that the probability of taking a setp to the left or to the
# right is equal. This yields a random direction of motion. However,
# Even if we mildly change the bias in the coin flip we can get something that
# is very reminiscent of directed motion. but let's prove it to ourselves!
# We'll start by defining some parameters.
n_steps = 1000
n_simulations = 100

# Set an array of step probabilities.
p = [0.5, 0.55, 0.6, 0.9]

# For this simulation, we will want to iterate through each probability,
# simulation, and step. This will have to be a three-layeted for loop.

# Make a storage array for the displacement.
displacement = np.zeros((len(p), n_simulations, n_steps))
for i in range(len(p)):
    for j in range(n_simulations):
        # We want to mark our initial position for each simulation.
        position = 0

        # Loop through each step and flip the coin.
        for k in range(n_steps):
            flip = np.random.rand()

            if flip < p[i]:
                position += 1
            else:
                position -= 1

            #  Now we'll store the displacement at each step.
            displacement[i, j, k] = position

# To generate the plot, we'll show the different walks on different axes on
# the same figure. This is something called a 'subplot'.
fig, ax = plt.subplots(2, 2)
ax = ax.ravel()  # Makes  it iterable.

# This will generate a figure with two rows of axes each with two columns. We
# can acccess these axes simply by iterating through the ax variable.

# We also want to set up a vector of steps so we know where each walker is.
time = np.arange(0, n_steps, 1)
for i in range(len(p)):
    for j in range(n_simulations):
        ax[i].plot(time, displacement[i, j, :])

    ax[i].set_title('p = ' + str(p[i]))
    ax[i].set_xlabel('time (number of steps)')
    ax[i].set_ylabel('position')
plt.show()
