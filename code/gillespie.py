# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this script, we will perform a simulation of constitutive expression from
# a promoter where the mRNA is allowed to decay. We are assuming that the mRNA
# is transcribed a constant rate. We'll begin by defining some parameters of
# our simulation.
prod_rate = 20  # Production rate of mRNA in units of 1/min
gamma = 1 / 3  # Degradation rate of mRNA in units of 1/min
dt = 0.01  # Time step in units of minutes.
total_time = 30  # Total length of the simulation in units of minutes.
time = np.arange(0, total_time, dt)  # Time vector

# We can do a simple euler integration of the mRNA transcription dynamics.
m_theo = np.zeros(len(time))  # A storage vector for the theoretical mean.
for i in range(1, len(time)):
    m_theo[i] = m_theo[i-1] + prod_rate * dt - gamma * dt * m_theo[i-1]

# For the gillespie simulation, we will determine whether an mRNA is
# transcribed based on the probabilities of production, decay, or no change at
# all at each point in time.
n_simulations = 100  # Total number of simulations to perform
max_it = 1000  # Maximum number of algorithm iterations.
m_sim = np.zeros((n_simulations, max_it))
tau = np.zeros((n_simulations, max_it))

for i in range(n_simulations):
    for j in range(1, max_it):
        degradation = gamma * m_sim[i, j-1]
        prob_sum = prod_rate + degradation

        # Determine the time to the next reaction (tau). We can pull this from
        # a Poisson distribution by relating it to a uniform distribution as
        # follows
        tau[i, j] = 1 / prob_sum * np.log(1 / np.random.rand())

        # Now we'll flip a coin and decide which step to take.
        flip = np.random.rand()
        if flip <= prod_rate / prob_sum:
            m_sim[i, j] = m_sim[i, j-1] + 1
        else:
            m_sim[i, j] = m_sim[i, j-1] - 1


# And our simulation is done! We need to now relate teh tau vector to physical
# time. We can do this simply by summing up each step of the tau.
plt.figure()

# Plot each simulation
for i in range(n_simulations):
    t = np.empty_like(tau[i])
    for j in range(len(tau[i])):
        t[j] = np.sum(tau[i, :j])

    # Now plot the simulation vs time.
    plt.plot(t, m_sim[i, :])

# Plot the thoeretical prediction.
plt.plot(time, m_theo, 'k-')

# Set the labels.
plt.xlabel('time (min)')
plt.ylabel('number of mRNAs')
plt.ylim([0, 100])
plt.show()
