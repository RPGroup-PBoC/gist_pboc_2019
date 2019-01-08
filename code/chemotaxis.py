# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this sscript, we will examine the probability of a chemoreceptor on a
# bacterium being active over a range of ligand concentrations. In the case of
# chemotaxis, the receptor is only active when the receptor is unbound. When
# bound with a particular ligand, the receptor becomes inactive, resulting in
# the bacterium swimmming forward. Methylation/demethylation of these receptors
# allows the bacterium to adapt to varying background levels of
# chemoattractant. THis chemical modification effectively changes the active
# energy of the receptor. Let's pretend that we are making those chemical
# modifications.

# We'll start by defining some parameters.
e_a = [-10, -4, -1]  # Range of active receptor binding energy (in kT).
e_i = 0  # Inactive receptor binding energy (in kT).
deb_a = -5  # Energetic difference of binding to an active receptor (in kT).
deb_i = -10  # Energetic difference of binding to an inactive receptor (in kT)
L = np.logspace(-5, -1, 1000)  # Range of ligand concentrations (in M).

# Since we are interested in a range of active energies, we'll evaluate our
# function at each e_a through a for-loop.
plt.figure()
for i in range(len(e_a)):
    # To save ourselves some typeing, we can write the partition function
    # explicityly.
    Z = np.exp(-e_a[i]) * (1 + L * np.exp(-deb_a))**2 +\
        np.exp(-e_i) * (1 + + L * np.exp(-deb_i))**2

    # Now that we ahve the partioion function, let's evaluate P_active.
    p_active = np.exp(-e_a[i]) * (1 + L * np.exp(-deb_a))**2 / Z

    # Now we can plot it. We'll add a legend at each step.
    plt.semilogx(L, p_active, label=r'$\epsilon_a$ = ' + str(e_a[i]) + ' kT')

# Most importantly...
plt.xlabel('ligand concentration (M)', fontsize=15)
plt.ylabel(r'$p_{active}$', fontsize=15)
plt.legend()
plt.show()

# The shifting of the curves demonstrates adaptation. Notice how the 'sweet
# spot' (inflection point) is always present, jsut at different points. This
# is how a bacterium can reverse direction over four orders of magnitude
# regarding ligand concentration.
