# of course...
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set_context('talk')
plt.close('all')

# Set up our parameters.
phi = np.linspace(0, 1, 500)
chi = [0.01, 0.1, 0.5, 0.8, 1, 1.25, 1.5, 2, 4]

# Make a figure, iterate over chi, and plot the free energy
plt.figure()
for x in chi:
    # Calculate the free energy
    free_energy = 2 * x * phi * (1 - phi) + phi * np.log(phi) + (1 - phi) * np.log(1 - phi)
    plt.plot(phi, free_energy, label='$\chi$ = ' + str(x))

# Add labels.
plt.xlabel('$\phi$')
plt.ylabel(r'$\frac{F(\phi)}{k_BT}$')
plt.title('Regular Solution Theory')
plt.tight_layout()
plt.savefig('../figures/regular_solution_theory.png', bbox_inches='tight')
plt.legend()
plt.show()

# Plot the Flory-Huggins theory.
L = 100 # Degree of polymerization
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
plt.tight_layout()
plt.savefig('../figures/flory_huggins_theory.png', bbox_inches='tight')
plt.show()
