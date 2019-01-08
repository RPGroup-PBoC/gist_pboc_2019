import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
plt.close('all')
# In this tutorial, we will investigate the dynamics of mRNA expression from an
# unregulated (constitutive) promoter.

# As we derived in class, the differential equation we wish to integrate is
#               dm/dt = r - y*m
# where r is the production rate, y is the degradation rate, and m is the
# number of mRNAs in the cell.
# Before we actually do the integration, we'll first examine how both the
# production term and the degradation term contribute to the change in mean
# mRNA copy number. To start, we will define some parameter values.

# Assign values to our model parameters.
r = 1  # mRNA production rate in 1/min.
gamma = 1/3  # mRNA decay rate in 1/min.
num_mRNA = np.arange(0, 10, 1)  # Set up a vector of mRNAs from 0 to 10

# Determine the contribution of each component.
prod = np.ones_like(num_mRNA) * r  # Vector of production component.
deg = gamma * num_mRNA  # Vector of degradation component.

# Plot the two together.
plt.figure()
plt.plot(num_mRNA, prod, '-', label='production')
plt.plot(num_mRNA, deg, '-', label='degradation')

# Add a legend and axis labels.
plt.legend(loc='upper left')
plt.xlabel('number of mRNA')
plt.ylabel(r'contribution to $dm/dt$')
plt.show()

# The point at which they intersect is the equilibrium steady state. We can
# solve our differential equation for this by setting the derivative equal to 0
# and solving for m. This yields
#                   dm/dt = 0 = r/y.
# Plugging in our values gives a value of 3, which is precisely at the
# intersection point of our plot. Now that we have a good sense at how each
# component contributes to the change in mRNA copy number, we can perform the
# integration.

# Define the parameters of our integration
time = 20  # in units of minutes.
dt = 0.1  # time step in units of min
num_steps = int(time / dt)  # the number of time steps to integrate.
time_vec = np.linspace(0, time, num_steps)  # Define the time vector.
init_condition = 0  # number of mRNA molecules at time 0.

# We will want to keep track of the mRNA at each time point, so we need to
# generate a storage vector.

# Set up a storage vector for our mRNA and set the initial condition.
m_t = np.empty_like(time_vec)
m_t[0] = init_condition

# Loop through each time point and evaluate the differential equation.
for i in range(1, num_steps):
    m_t[i] = m_t[i-1] + r*dt - gamma * m_t[i-1]*dt

# Plot the result of our integration.
plt.figure()
plt.plot(time_vec, m_t, '-', label='integration')

# Plot the steady state as a horizontal line.
plt.hlines(r/gamma, xmin=0, xmax=time, label='steady state')

# Extend the margins and add axis labels.
plt.ylim([0, 5])
plt.xlabel('time (min)')
plt.ylabel(r'$m(t)$')
plt.legend()
plt.show()

# We can see that after about 15-ish minutes, the system has reached the
# equilibrium steady state. What happens when we change our initial condition?
# Let's consider the case where there are 10 mRNAs to begin with in the cell.

# Set up a storage vector for our mRNA and set the new initial condition.
m_t2 = np.empty_like(time_vec)
m_t2[0] = 10  # Our new initial condition

# Loop through each time point and evaluate the differential equation.
for i in range(1, len(time_vec)):
    m_t2[i] = m_t2[i-1] + r*dt - gamma * m_t2[i-1]*dt

# Now we'll plot everything on the same axes.
plt.figure()
plt.plot(time_vec, m_t, '-', label=r'$m_0 = 0$')
plt.plot(time_vec, m_t2, '-', label=r'$m_0 = 10$')
plt.hlines(r/gamma, xmin=0, xmax=time, label='steady state')

# Change the axes limits and add labels.
plt.ylim([0, 12])
plt.xlabel('time (min)')
plt.ylabel(r'$m(t)$')
plt.legend()
plt.show()

# In this tutorial, we have learned a bit about mRNA copy number dynamics for a
# case of simple constitutive expression. With a few modifications, we can
# include some simple modes of regulation while making a few changes to our
# differential equation. Additionally, we learned the very useful technique of
# numerical integration using the Euler method.
