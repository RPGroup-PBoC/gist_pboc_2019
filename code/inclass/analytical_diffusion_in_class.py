# For doing math
import numpy as np

# For plotting
import matplotlib.pyplot as plt
import seaborn

# Define our diffusion constant.
D = 10  # in microns squared per second

# Time for diffusion to occur
time_steps = [1, 5, 10, 50]

# Set up an array of distances
x = np.linspace(-100, 100, 1000)

# Plot the probabilities
plt.figure()
for t in time_steps:
    # Compute and plot!
    p_t = 1 / np.sqrt(4 * np.pi * D * t) * np.exp(-x**2 / (4 * D * t))

    # plot
    plt.plot(x, p_t, label=t)

plt.xlabel('distance in microns')
plt.ylabel('P_t')
plt.legend()
plt.show()
