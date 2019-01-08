# Import our favorite modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# Define some parameter values.
r = 1  # mRNA production rate in 1/min
gamma = 1 / 3  # mRNA decay in 1/min
time = 20  # in min
dt = 0.1  # time step in min
num_steps = time / dt
init_cond = 10  # in units of number of mRNA
time_vec = np.linspace(0, time, num_steps)

# Make a vector to store the mRNA count
m_t = np.zeros_like(time_vec)
m_t[0] = init_cond
# Do the integration!
for t in range(1, int(num_steps)):
    m_t[t] = m_t[t-1] + r * dt - gamma * dt * m_t[t-1]

plt.figure()
plt.plot(time_vec, m_t, 'r-', label='m(0) = ' + str(init_cond))
plt.xlabel('time (min)')
plt.ylabel('$m(t)$')
plt.ylim([0, 10])

# Do the integration again with a different initial condition.
m_t = np.zeros_like(time_vec)
m_t[0] = 0
# Do the integration!
for t in range(1, int(num_steps)):
    m_t[t] = m_t[t-1] + r * dt - gamma * dt * m_t[t-1]


plt.plot(time_vec, m_t, 'b-', label='m(0) = 0')
plt.legend()
plt.show()
