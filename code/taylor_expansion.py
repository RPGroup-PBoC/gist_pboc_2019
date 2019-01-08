# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# This script will demonstrate the usefulness of using the taylor expansion
# to approximate the local shape of a curve. We will look at a simple function,
# y = cos(x), but it may be useful for you to try more complicated functions.

# To start, we'll make an array of x values at which we will evaluate cos(x).
x = np.linspace(-4 * np.pi, 4 * np.pi, 500)
y = np.cos(x)

# Let's evaluat eh first handful of orders.
order_0  = np.ones_like(x)
order_2 = 1 - x**2 / np.math.factorial(2)
order_4 = order_2 + x**4 / np.math.factorial(4)
order_6 = order_4 - x**6 / np.math.factorial(6)

# Now let's plot all of them together.
plt.figure()
plt.plot(x, y, label='cos(x)')
plt.plot(x, order_0, label='0th order')
plt.plot(x, order_2, label='2nd order')
plt.plot(x, order_4, label='4th order')
plt.plot(x, order_6, label='6th order')

# Add the legend and labels.
plt.legend()
plt.xlabel('x')
plt.ylabel('y')

# Change the yrange so it's easier to see.
plt.ylim([-1.5, 1])
plt.show()

# With this, we can see that even the first few expansions are a great
# approximation of the actual behavior.
