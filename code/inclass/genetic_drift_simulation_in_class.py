# IMPORT OUR THINGS
import numpy as np
import matplotlib.pyplot as plt
import seaborn

plt.close('all')
# Set up some parameters.
N = 10
num_gen = 200

# Our initial allele frequencies.
p = 0.5  # Our initial 'A' frequency

# Empty storage list.
freq_a = []
for i in range(num_gen):
    # Flip our coins.
    flip = np.random.rand(N)

    # Count the number of A alleles.
    num_a = np.sum(flip < p)

    # Compute the frequency.
    freq = num_a / len(flip)

    # Store the frequency.
    freq_a.append(freq)
    p = freq

# Make a vector of generations.
gen_vec = np.arange(0, num_gen, 1)

# Plot it!
plt.figure()
plt.plot(gen_vec, freq_a)
plt.xlabel('number of generations')
plt.ylabel('frequency of A')
plt.show()
