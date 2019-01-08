# Import the necessary things.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pboc_utils as pboc

import scipy.misc
# Define some parameters
N = 8  # number of individuals in the population
P = np.zeros((2 * N + 1, 2 * N + 1))

# Populate the matrix
for i in range(2 * N + 1):
    for j in range(2 * N + 1):
        P[j, i] = scipy.misc.comb(2 * N, j) * (i / (2 * N))**j *\
            (1 - i / (2 * N))**(2 * N - j)


# Now apply it to the generation vector
n_gen = 80
allele_freq = np.zeros((2 * N + 1, n_gen))
allele_freq[N, 0] = 1.0

# Apply the matrix on each generation
for i in range(1, n_gen):
    allele_freq[:, i] = np.dot(P, allele_freq[:, i - 1])

pboc.bar3(allele_freq, xlabel='number of generations', ylabel='number of alleles', zlabel='probability')
plt.show()
