# import the important modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io

# This script will simulate diffusion of a particle over time. While we can
# solve the diffusion equation analytically, this will serve as a fun way to
# learn how to write stochastic simulations.

# We'll begin by considering the diffusion of a single particle in one
# dimension. We will say that at any point in time, the particle has an equal
# probability of taking a step to the left as it can take a step to the right.
# Let's takea  look at the trajectory of a single particle.

# We'll seed the random number generator such that our results will be
# the same each time we run it.
#np.random.seed(42)

# Define the nubmer of steps the simulation will run.
n_steps = 1000

# Define the stepping probability.
p = 0.5  # Equal probability of stepping right or left.

# We can simulate the random walk by flipping a coin at each step and testing
# its value with the probability.

position = 0
displacement = np.zeros(n_steps)

for i in range(n_steps):
    # Flip a coin.
    flip = np.random.rand()

    # Test what the value of this flip is.
    if flip < p:
        position = position + 1
    else:
        position = position - 1

    # We want some way to store the position of our walker at this step.
    displacement[i] = position


# And that's it! Let's plot eh trajectory as a function of time (step number)
# and see where it went.

plt.figure()
time = np.arange(0, n_steps, 1)
plt.plot(time, displacement)
plt.xlabel('time (number of steps)')
plt.ylabel('position')
plt.show()

# It's not very useful to see where a single particle diffused. Let's up our
# number of walkers to something huge. That will be representative of the
# behavior of our system. Let's modify our above loop to do this.
n_simulations = 1000
displacement = np.zeros((n_simulations, n_steps))
for i in range(n_simulations):

    # Reset the position for each simulation
    position = 0

    # Loop through each step.
    for j in range(n_steps):

        # Flip a coin.
        flip = np.random.rand()
        if flip < p:
            position = position + 1
        else:
            position = position - 1

        # Keep track of the displacement for each walker.
        displacement[i, j] = position

# Now, let's loop through each simulation and plot the trajectory.
plt.figure()
for i in range(n_simulations):
    plt.plot(time, displacement[i, :])

plt.xlabel('time (number of steps)')
plt.ylabel('position')
plt.show()

# Now we can see a bunch of different trajectories. Let's look at the
# distribution of ending positions. We can make a new array using the
# sotred displacements from the previous simulation.
ending_position = np.zeros(n_simulations)
for i in range(n_simulations):
    ending_position[i] = displacement[i, -1]


# Plot the distribution.
plt.figure()
plt.hist(ending_position, bins=100)
plt.xlabel('ending position')
plt.ylabel('counts')
plt.show()


# Let's see how the displacemtn would change as a function of the step number.
# In our simulation, we are storing positions that are both positive and
# netative. If we wanted to get an idea of the average distance traveled, we
# can simply square the final position and take the average of all of the
# simulations.
n_steps = [10, 100, 500, 1000, 5000, 10000]
final_position = np.zeros(n_simulations)
msd = np.zeros_like(n_steps)
for i in range(len(n_steps)):
    for j in range(n_simulations):
        position = 0
        for k in range(n_steps[i]):
            # Perorm the coin flip.
            flip = np.random.rand()
            if flip < p:
                position = position + 1
            else:
                position = position - 1

        # Save the final position.
        final_position[j] = position

    # Compute and store the mean square displacement.
    msd[i] = np.mean(final_position**2)


# Let's look at how it scales with time.
plt.figure()
plt.loglog(n_steps, np.sqrt(msd), 'o', label='simulation')

# We predict the root mean squared displacement would scale linear with the
# root of the number of steps we've taken
n_vec = np.logspace(-1, 5, 500)
plt.loglog(n_vec, np.sqrt(n_vec), '-', label='data')
plt.xlabel('time (number of steps)')
plt.ylabel('mean squared displacement')
plt.show()
