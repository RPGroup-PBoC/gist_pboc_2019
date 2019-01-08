# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# In this script, we will examine the behavior fo the fold-change equation for
# the case of the LacI inducible regulation. As a reminder, the equation for
# fold-change is given by
#
# fold-change = (1 + (1 + c/ka)^2 / ((1 + c/ka)^2 + e^(-ep_ai)(1 + c/ki)^2) * (R / Nns)e^(-ep_r))^-1
#
# where c is the inducer concentration, ka is the binding constant of inducer
# to the active repressor, ki is the binding constant of inducer to the
# inactive repressor, ep_ai is the energy difference between active and
# inactive repressor, R is the number of repressors in the cell, Nns is the
# number of nonspecific binding sites for the repressor and ep_r is the
# binding energy of the repressor to the DNA.

# We'll look at how this changes under three circumstances -- when there is
# varying numbers of repressors per cell, with varying binding strengths of the
# repressor to the DNA, and with varying energetic difference between the
# active and inactive state.

# To begin, we'll define some parameter values for varying the repressor copy
# number.
c = np.logspace(-9, -2, 500)  # Range of inducer concentration in M
ka = 200E-6  # Binding constant of inducer to repressor in M
ki = 50E-9  # Binding constant of inducer to inactive repressor in M
ep_ai = -5  # Difference between active and inactive state in kT
ep_r = -14  # Binding energy of the repressor to the DNA in the active state in kT

# Now we'll set up a range of R.
repressors = [10, 50, 100, 500, 1000]

# We'll use the same equation for fold-change over and over again, so we
# can define this as a function.
def calc_fold_change(c, ka, ki, ep_ai, ep_r, R, Nns=5E6):
    """
    Computes the fold change of the MWC model for lacI induction.
    """

    # Compute the mwc component
    mwc_term = (1 + c / ka)**2 / ((1 + c/ka)**2 +\
                np.exp(-ep_ai) * (1 + c /  ki)**2)

    # Compute the repression.
    repression = 1 + mwc_term * (R / Nns) * np.exp(-ep_r)

    # Now compute and return the fold-change
    fold_change = 1 / repression
    return fold_change


# Great! Now let's evaluate this function over a range of repressor copy
# numbers.
plt.figure()
for R in repressors:
    fold_change = calc_fold_change(c, ka, ki, ep_ai, ep_r, R)
    plt.plot(c, fold_change, label='R = ' + str(R))

# Now we can add our labels.
plt.xscale('log')
plt.xlabel('inducer concentration (M)')
plt.ylabel('fold-change')
plt.legend()
plt.show()

# We can see that with a few repressors, it takes only a little bit of inducer
# to greatly increase the expression.

# Let's take a look at a range of DNA binding strengths.
R = 250  # Reset the number of repressors as static.
binding_strengths = [-5, -7, -10, -15, -18, -20]  # Range of DNA binding strengths in kT

# Set up a new figure and loop once again.
plt.figure()
for ep_r in binding_strengths:
    fold_change = calc_fold_change(c, ka, ki, ep_ai, ep_r, R)
    plt.plot(c, fold_change, label=r'$\Delta\epsilon_R$ = ' + str(ep_r) + ' kT')

# Add labels and show it.
plt.xscale('log')
plt.xlabel('inducer concentration (M)')
plt.ylabel('fold-change')
plt.legend()
plt.show()

# Now we can see that by changing that binding sterength, we can see a large
# amount of diversity of response. We can significantly change the dynamic
# range of fold-change.

# Now let's take a look at varying the energy difference between the active and
# inactive repressor.
ai_energies = [-8, -5, -3, -2, -1, 0]
ep_r = -14

# Set up a new figure and loop once again... again.
plt.figure()
for ep_ai in ai_energies:
    fold_change = calc_fold_change(c, ka, ki, ep_ai, ep_r, R)
    plt.plot(c, fold_change, label=r'$\Delta\epsilon_{ai}$ = ' + str(ep_ai) + ' kT')
plt.xscale('log')
plt.xlabel('inducer concentration')
plt.ylabel('fold-change')
plt.legend()
plt.show()
