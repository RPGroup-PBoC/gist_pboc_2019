# For doing math
import numpy as np

# For plotting
import matplotlib.pyplot as plt
import seaborn

# Define our diffusion constant.
D = 10  # in square microns per second

# Generate an array of lengths.
L = np.logspace(-3, 6, 500)

# Determine the diffusion time.
time = L**2 / (2 * D)

plt.figure()
plt.loglog(L, time)
plt.xlabel('length in microns')
plt.ylabel('time in seconds')
plt.show()
