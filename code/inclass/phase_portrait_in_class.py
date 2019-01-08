# Duhhhh
import numpy as np
import matplotlib.pyplot as plt
import seaborn

plt.close('all')

# Define the parameters
r = 20  # the production rate
gamma = 1 / 30  # the degradation rate
k = 200 # in units of concentration
max_R = 1000  # maximum number of R1 and R2
R1 = np.linspace(0, max_R, 500)
R2 = np.linspace(0, max_R, 500)

# Compute the nullclines.
R1_null = (r / gamma) / (1 + (R2 / k)**2)
R2_null = (r / gamma) / (1 + (R1 / k)**2)

# Plot the nullclines.
plt.figure()
plt.plot(R1, R1_null, label='dR1/dt = 0')
plt.plot(R2_null, R2, label='dR2/dt = 0')
plt.xlabel('R1')
plt.ylabel('R2')
plt.legend()
plt.show()

# Generate the vector fields
R1_m, R2_m = np.meshgrid(R1[1::30], R2[1::30])

# Compute the derivatives
dR1_dt = -gamma * R1_m + r / (1 + (R2_m / k)**2)
dR2_dt = -gamma * R2_m + r / (1 + (R1_m / k)**2)

# Plot the vector fields!!
plt.quiver(R1_m, R2_m, dR1_dt, dR2_dt)
plt.show()


# Plot the orbit.
time = 200
R1 = 800
R2 = 400

# Loop through time and integrate.
for t in range(time):
    dR1 = -gamma * R1 + r / (1 + (R2 / k)**2)
    dR2 = -gamma * R2 + r / (1 + (R1 / k)**2)

    # Add this change to our current position
    R1 = R1 + dR1
    # This is the same operation as above..
    R2 += dR2

    plt.plot(R1, R2, 'ro')
    plt.show()
    plt.pause(0.05)
