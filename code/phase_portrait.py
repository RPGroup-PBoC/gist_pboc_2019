# Import the necessary modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# In this script, we will generate the null clines and the phase portrait
# for our genetic switch circuit. As a reminder, consider we have two molecules
# (R1 and R2) which repress each other. We can summarize the dynamics of this
# system by writing the differential equations of
#
#   dR1/dt = -gamma * R1 + (r / (1 + (R1 / k)^2))
#   dR2/dt = -gamma * R2 + (r / (1 + (R2 / k)^2))
#
# where gamma is the degradation rate of the proteins, r is the production
# rate of the proteins, and k is the dissociation constant of the repressor
# for the DNA. In this script, we are assuming that the production rate, the
# degradation rate, and the dissociation constant are the same for both R1 and
# R2. To begin, we will define some parameter values.

r = 20  # production rate of the proteins in units of 1 / seconds
gamma = 1 / 30 # degradation rate of proteins in units of / seconds
k = 200 # dissociation constant of the repressor for the DNA in M.
num_R = 800
R1 = np.linspace(0, num_R, 500)
R2 = np.linspace(0, num_R, 500)

# It will be useful to look at the nullclines of this system. The nullclines
# are defined as the points in which the derivatives across the two components
# of our system is equal to zero. We can solve these equations written above
# and generate the following solutions:
R1_null = (r / gamma) * 1 / (1 + (R2 / k)**2)
R2_null = (r / gamma) * 1 / (1 + (R1 / k)**2)

# With the nullclines defined, let's plot them and give them proper labels.
plt.figure()
plt.plot(R1, R1_null, label=r'$\frac{dR1}{d}t$ = 0')
plt.plot(R2_null, R2, label=r'$\frac{dR2}{dt}$ = 0')
plt.xlabel('R1')
plt.ylabel('R2')
plt.legend()


# We can see that there are three points at which the nullclines intersect.
# these are 'fixed points'. The two points at the upper left and lower right
# are 'stable' and the one right in the middle is unstable. We can prove this
# to ourselves by evaluating the derivates at a large array of points and
# plotting the vectors. For this, we will need to generate a matrix of values
# to compute the derivatives over. We can do this through the numpy meshgrid
# function.
R1_m, R2_m = np.meshgrid(R1[::30], R2[::30])  # only generating every 30th point

# Now we'll evaluate the derivatives on this mesh.
dR1_dt = -gamma * R1_m + r / (1 + (R2_m / k)**2)
dR2_dt = -gamma * R2_m + r / (1 + (R1_m / k)**2)

# Now plot the vector fields.
plt.quiver(R1_m, R2_m, dR1_dt, dR2_dt)
plt.show()

# We can see that the vectors point towards the stable fixed point and always
# from the unstable point. If we choose a few random positions in the parameter
# space, we should be able to 'watch' the integration take place. We'll
# randomly choose 10 positions to start from.

# Choose a random starting position.
N = 10
R1_init = np.random.rand(N) * num_R
R2_init = np.random.rand(N) * num_R

# We'll integrate over a total of 200 time steps.
total_time = 200

# We'll set the initial condition.
R1 = R1_init
R2 = R2_init
dt = 1 # Time step. This is 1/30th of 1/gamma

# Now we'll iterate over time and compute the change.
for t in range(total_time):
    dR1 = -gamma * dt * R1 + r * dt / (1 + (R2 / k)**2) * dt
    dR2 = -gamma * dt * R2 + r * dt / (1 + (R1 / k)**2) * dt

    # Now add the change to our current position.
    R1 += dR1
    R2 += dR2

    # Now we'll plot the curve.
    plt.plot(R1, R2, 'r.')
    plt.show()

    # Now we'll tell it to wait
    plt.pause(0.005)

# Cool! We can see the integration take place and watch points migrate away
# from the unstable fixed point and towards the stable one.s
