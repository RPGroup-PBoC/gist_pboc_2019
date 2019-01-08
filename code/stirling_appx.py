# IMport the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this short script, we will see how Stirling's approximation can be
# useful for approximating very large numbers. For example, computing 200! is
# very difficult because that is beyond the basic precision of a 64-bit
# computer. We can do better by converting the numbers to something called
# 'arbitrary precision', but it is often sufficient to perform approximations
# by hand.

# As a reminder, Stirling's approximation is
#
#     log(n!) = n * log(n) - n
#
# We'll show this by computing the log of the factorial for large numbers
# and show it along with the approximation.

# We'll first set upa range of n to compute.
n_array = np.arange(0, 150, 1)

# Now we'll compute the factorial and store the log.
log_factorial = []
for n in n_array:
    n_factorial = float(np.math.factorial(n))   # Converts to a float for precision
    log_factorial.append(np.log(n_factorial))

# Now we can compute our approximation.
stirlings_approx = n_array * np.log(n_array) - n_array

# Now let's plot the log factorial and our approximation together.
plt.figure()
plt.plot(n_array, log_factorial, label='factorial')
plt.plot(n_array, stirlings_approx, label='approximation')
plt.xlabel('$n$')
plt.ylabel('$n!$')
plt.legend()
plt.show()

# We see that as n becomes larger and larger, our approximation gets even
# better.
