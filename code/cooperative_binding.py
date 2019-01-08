# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this script, we will consider a receptor which can bind two ligands.
# Using the statistical mechanics in we learned in lecture, we can enumerate
# these states, figure out the multiplicity, and compute their statistical
# weights.  How does the probability of each state change with a concentration
# of ligand? Let's find out!

# We'll start by declaring some variables.
d_eb = -5  # Energetic difference between bound and unbound states (in kT)
d_ei = -2  # Energetic difference of interaction
L = np.logspace(-4, -1, 1000)  # Range of ligand concentration.


# To save our fingers some work, we'll define the partition function.
Z = 1 + 2 * L * np.exp(-d_eb) + L**2 * np.exp(-(2 * d_eb + d_ei))

# Now we'll compute the probability of each state
p_empty = 1 / Z  # Nothing bound -- empty receptor.
p_single = L * np.exp(-d_eb) / Z  # Probability of one ligand being bound.
p_double = L**2 * np.exp(-(2 * d_eb + d_ei)) / Z  # Probability of both bound.

# Now we'll plot the probabilites as a function of ligand concentration.
plt.figure()
plt.semilogx(L, p_empty, '-', label='empty')
plt.semilogx(L, p_single, '-', label='singly bound')
plt.semilogx(L, p_double, '-', label='doubly bound')

# Add labels as one always should.
plt.xlabel('ligand concentration (M)')
plt.ylabel('probability of state')
plt.legend()
plt.show()


# We can see that with a high concentration of ligand, the most common state
# is having both bound at once. In the middle of the range, empty and doubly
# bound are equal with a singly bound state less probable. Let's now compute
# the average number of bound ligands as a function of concentration.

# Since we have already enumerated our probabilities, the calculation is
# trivial.
n_bound = 2 * p_single + 2 * p_double

# Plot it!
plt.figure()
plt.plot(L, n_bound, 'k-')
plt.xlabel('ligand concentration (M)')
plt.ylabel('average number of bound ligands')
plt.show()

# Now go back and try changing the energy of interaction and see how the
# curves change!
