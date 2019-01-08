# Import our favorite modules
import numpy as np
import matplotlib.pyplot as plt
# import seaborn
# Define our parameters
Nns = 5E6  # in units of number of binding sites
R = np.logspace(0, 3, 500)
de_r = -15  # in kT

# Calculate the fold_change
fold_change = 1 / (1 + (R / Nns) * np.exp(-de_r))

plt.figure()
plt.plot(R, fold_change)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('R')
plt.ylabel('fold-change')
plt.show()

# Define more parameters
P = 1E3  # Number of polymerases
de_p = -8  # in units of kT

# Write the denominator
denom = 1 + (R / Nns) * np.exp(-de_r) + (P / Nns) * np.exp(-de_p)

# Compute the probabilities
p_empty = 1 / denom
p_repressor = (R / Nns) * np.exp(-de_r) / denom
p_polymerase = (P / Nns) * np.exp(-de_p) / denom

plt.figure()
plt.plot(R, p_empty, label='empty')
plt.plot(R, p_repressor, label='repressor bound')
plt.plot(R, p_polymerase, label='polymerase bound')
plt.xscale('log')
# plt.yscale('log')
plt.xlabel('number of repressors')
plt.ylabel('probability of state')
plt.title('viral promoter binding energy')
plt.legend()
plt.show()
