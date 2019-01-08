# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this script, we will explore the dynamics of phase transitions in a
# two component mixture. We discussed two models in class -- the Regular
# Solution theory and the Flory-Huggins Theory.

# The regulular solution describes a mixture of componnents in which each
# individual molecule occupies only one lattice site in our closed system.
# We can define the probability of observing like-bonded molecules as phi.
# We can also include the Flory parameter which describes the energy of mixing
# the two species. The full description of the regular solution theory
# description of the free energy is
#
# free-energy / kBT = (2 * X * phi) * (1 - phi) + phi * log(phi) * (1-phi) * log(phi)#
#
# where X is the flory parameter. Let's see how varying this parameter
# influences the phases of the system.

# We'll set up our parameters.
phi = np.linspace(0, 1, 500)
chi = [0.01, 0.1, 0.5, 0.8, 1, 1.25, 1.5, 2, 4]

# We'll compute the value at each chi and plot it at once.
plt.figure()
for x in chi:
    # Calculate the free energy
    free_energy = 2 * x * phi * (1 - phi) + phi * np.log(phi) + (1 - phi) * np.log(1 - phi)
    plt.plot(phi, free_energy, label='$\chi$ = ' + str(x))

# Of course we should add labels.
plt.xlabel('$\phi$')
plt.ylabel(r'$\frac{F(\phi)}{k_BT}$')
plt.title('Regular Solution Theory')
plt.legend()
plt.show()

# We see that at the intermediate values of X , there are now multiple local
# minima where phases will separate.


# Now we can consider a case in which we are mixing together polymers. In this
# case, the polymers now occupy several different lattice sites. The result for
# the free energy of the system is given by the Flory-Huggins theory,
#
# free-energy / k_BT = 2 * X * phi * (1-phi) + (phi / L) * log(phi) + (1-phi) * log(1-phi)
#
# where L is the degree of polymerization, essentially just the length of the
# polymer.
L = 100 # Degree of polymerization

# Now we'll do the same iteration of x and  plot the free energy.
plt.figure()
for x in chi:
    # Calculate the free energy.
    free_energy = 2 * x * phi * (1 - phi) + (phi / L) * np.log(phi) +\
            (1 - phi) * np.log(1 - phi)
    # Plot it!
    plt.plot(phi, free_energy, label='$\chi$ = ' + str(x))

plt.xlabel('$\phi$')
plt.ylabel(r'$\frac{F(\phi)}{k_BT}$')
plt.legend()
plt.title('Flory-Huggins Theory')
plt.show()

# We now see that there is a degree of nonsymmetry when compared to the
# regular solution theory. This significantly changes the way that
# phases separate.
