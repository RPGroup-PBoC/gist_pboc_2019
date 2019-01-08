# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# In this script we will look at the time to diffuse a given distance as
# a function of time. To begin, let's define our diffusion constant.
D = 10**7  # in units of micron squared per second.

# Now we'll set up an array of distances.
L = np.logspace(0, 8, 1000)

# We know that the time to diffuse is related to the length as
#
#       time = (L**2) / (2 * D)
#
#  where L is the length and D is the diffusion constant.
time = L**2 / (2 * D)

# Let's plot it!
plt.figure()
plt.loglog(L, time)
plt.xlabel('distance (in microns)')
plt.ylabel('time (in seconds)')
plt.show()

# From this plot, we can see that it would take a little over a year to
# diffuse a single centimeter!
