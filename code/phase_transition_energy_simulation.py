# Import our necessary things
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.close('all')
# As an example, let's generate our matrix.
size = 8
lattice = np.random.rand(size, size)  # Full of random numbers between 0 and 1

# Set our proportion of A and B.
p_a = 0.1
p_b = 1 - p_a

# Populate the matrix
pop_lattice = lattice < p_a

# Compute the energy.
plt.matshow(pop_lattice)
plt.show()

# Loop over to find the number of AA bonds.
num_aa = 0
for i in range(size - 2):
    for j in range(size - 2):
        sum_z = np.sum(pop_lattice[i:i+2, j:j+2])
        if sum_z > 1:
            num_aa += sum_z

print(num_aa)
