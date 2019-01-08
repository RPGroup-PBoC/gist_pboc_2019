# Import our modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# Define our parameter values.
c = np.logspace(-9, -2, 500)
ka = 150E-6   # Binding constant for inducer to active repressor in M
ki = 50E-9   # Binding constant for inducer to inactive repressor in M
ep_r = -14  # binding energy of repressor to DNA in units of kT
repressors = [10, 50, 100, 1000]  # Range of repressors
Nns = 5E6  # number of nonspecific binding sites
ep_ai = -5 # difference between active and inactive state in kT
# Loop and plot.
plt.figure()
for R in repressors:
    mwc_term = (1 + c / ka)**2 / ((1 + c/ka)**2 + np.exp(-ep_ai) * (1 + c/ki)**2)
    repression = 1 + mwc_term * (R / Nns) * np.exp(-ep_r)
    fold_change = 1 / repression

    # Plot and label it.
    plt.semilogx(c, fold_change, label='R = ' + str(R))

plt.xlabel('inducer concentration (M)')
plt.ylabel('fold-change')
plt.legend()
plt.show()
