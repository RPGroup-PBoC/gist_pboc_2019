# Import our stuff
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set_context('talk')

# Set up an array of n values
n_array = np.arange(0, 150, 1)

# Compute the factorial.
n_factorial = []
for n in n_array:
    value = float(np.math.factorial(n))
    n_factorial.append(value)

# Compute the log of the factorial
log_factorial = np.log(n_factorial)

# Compute Stirling's approximation.
stirling = n_array * np.log(n_array) - n_array

# Make the plot.
plt.figure()
plt.plot(n_array, log_factorial, label='log(n!)')
plt.plot(n_array, stirling, label='nlog(n) - n')
plt.xlabel('n')
plt.ylabel('log(n!)')
plt.legend(loc='upper left')
plt.tight_layout()
plt.savefig('../figures/stirlings_approximation.png', bbox_inches='tight')
plt.show()
