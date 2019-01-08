# code to compute cell growth


# Load the plotting utility matplotlib with an alias plt.
#import numpy as np
#import matplotlib.pyplot as plt



# Set the initial number of cells in the experiment.
number_of_cells = 1

# Set the number of division cycles for the experiment.
number_of_divisions = 10

# Set a list of number of cells at division d and start with 1 cell
N_d = [number_of_cells]
# Loop through each division event
for i in range(number_of_divisions):
    # Make the cells duplicate
    number_of_cells = number_of_cells * 2
    # Add the new number of cells to our storage list.
    N_d.append(number_of_cells)

# Print the result of our simulation.
print(N_d)

# Establish a vector of the division cycles
division_vector = np.arange(0, number_of_divisions + 1, 1)
print(division_vector)

# Generate the plot.
plt.plot(division_vector, N_d, 'o', label='simulation')
# Set the axis labels.
plt.xlabel('number of divisions')
plt.ylabel('number of cells')
# Add a legend
plt.legend()
