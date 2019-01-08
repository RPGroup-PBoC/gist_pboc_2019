# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Define the parameters
L = np.logspace(-12, 2, 1000)  # Range of ligand concentration.
k_da = 1E-9  # Binding constant for ligand to the active conformation in M
k_di = 100E-9   # Binding constant for ligand to the inactive conformation in M
ep_a = 5  # Energy of the active receptor in kT
ep_i = 0  # Energy of the inactive receptor in kT

# To save some typing, let's write out the partition function.
Z = np.exp(-ep_a) * (1 + L / k_da) + np.exp(-ep_i) * (1 + L / k_di)

# Now we can write the probabilites of each state.
p_a = np.exp(-ep_a) / Z  # Probability of the active unbound state
p_i = np.exp(-ep_i) / Z  # Probability of the inactive unbound state
p_ab = np.exp(-ep_a) * (L / k_da) / Z  # active bound state
p_ib = np.exp(-ep_i) * (L / k_di) / Z  # active bound state


# Now let's plot the probabilities.
plt.figure()
plt.plot(L, p_a, label='active unbound')
plt.plot(L, p_i, label='inactive unbound')
plt.plot(L, p_ab, label='active bound')
plt.plot(L, p_ib, label='inactive bound')
plt.xscale('log')
# plt.yscale('log')
plt.xlabel('ligand concentration (M)')
plt.ylabel('probability of state')
plt.legend()
plt.show()

# We can also plot the probability of a ligand being bound.
p_La = (L / k_da) / (1 + (L / k_da) + (L / k_di))
p_Li = (L / k_di) / (1 + (L / k_da) + (L / k_di))
plt.figure()
plt.plot(L, p_La, label='active bound')
plt.plot(L, p_Li, label='inactive bound')
plt.xlabel('ligand concentration (M)')
plt.ylabel('$P_{bound} (L)$')
plt.xscale('log')
plt.legend()
plt.show()
