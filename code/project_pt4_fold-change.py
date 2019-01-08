# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context('talk')
# For image processing.
import skimage.io

# For file management
import glob

# Utilities for our course.
import pboc_utils as pboc


# In the final part of our project, we will use the image processing tools
# we've developed to iterate over many different images and quantitatively
# compute the fold change in gene expression for a variety of different
# repressor copy numbers and operator binding energies.

# As a reminder, we are looking to test our statistical mechanical theory of
# gene expression. By enumerating all of the microscopic states and the
# statistical weights (given by the Boltzmann distribution), we arrive at
# the prediction,
#
#       fold-change = (1 + (R/Nns) * exp(-dEp / kBT))^-1,
#
# where R is the repressor copy number, Nns is the number of nonspecific
# binding sites available to the lacI repressor, dEp is the difference in
# energy between the DNA bound and unbound states of the repressor, and
# kBT is the thermal energy of the system. To be true biophysicists, let's plot
# our prediction first.

# Define some parameters of the prediction
ops = ['O1', 'O2', 'O3']
binding_energy = [-15.3, -13.9, -9.3]  # in units of kBT
R = np.logspace(0, 4, 500)
Nns = 4.6E6


# Define a function computing the fold change.
def fold_change(R, binding_energy, Nns):
    """
    Computes fold change.
    """
    return 1 / (1 + (R/Nns) * np.exp(-binding_energy))


# Generate a figure and plot the prediction for all operators.
plt.figure()

# Since we'll eventually plot data over these predictions, we'll make sure
# they will all be the correct colors.
colors = ['r', 'g', 'b']

# We'll loop through each binding energy and plot the result.
for i in range(len(binding_energy)):
    fc = fold_change(R, binding_energy[i], Nns)
    plt.plot(R, fc, '-', color=colors[i], label=ops[i] + ' theory')


plt.xlabel('number of repressors')
plt.ylabel('fold-change')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()
plt.tight_layout()
plt.savefig('figures/other_predictions.png', bbox_inches='tight')
# We see that we have lower fold-change (higher repression) when we have more
# repressor around to do the job. We also get lower fold-change if we increase
# the strength at which the repressor molecule binds the DNA. Let's put those
# image processing skills to work and test our theory. Experimentally, we
# will compute the fold-change as
#
#       fold-change = (<I_R> - <I_auto>) / (<I_delta> - <I_auto>)
#
# where <I_R> is the mean intensity in the presence of repressor, <I_auto> is
# the mean intensity of autofluorescence, and <I_delta> is the mean intensity
# in the absence of any repressor. Let's define the important part of the
# file names and the repressor copy numbers for each file first.
rep_names = ['auto', 'delta', 'wt', 'R60', 'R124', 'R260', 'R1220', 'R1740']
rep_numbers = [11, 60, 124, 260, 1220, 1740]  # This starts at 'wt' strain

# Note that the autofluorescence strain and the delta fluorescence strain are
# the same for each individual operator. To make things easier on ourselves,
# let's calculate those values first.

# Set up empty storage vectors for the means.
auto_means = []
delta_means = []

# Define the thresold value for our segmentation.
thresh = -0.2

# Loop through each operator.
for i in range(len(ops)):

    # Loop through the auto and delta strains. They are the first two in our
    # rep_names list.
    for j in range(2):

        # Get the names of all of the phase contrast files.
        phase_names = glob.glob('data/lacI_titration/' + ops[i] + '_' +
                                rep_names[j] + '_phase*.tif')
        fluo_names = glob.glob('data/lacI_titration/' + ops[i] + '_' +
                               rep_names[j] + '_yfp*.tif')

        # Set up a storage vector to store all of the single-cell intensities.
        intensities = []

        # Iterate through each pair of images and perform the image processing.
        for k in range(len(phase_names)):
            phase_im = skimage.io.imread(phase_names[k])
            fluo_im = skimage.io.imread(fluo_names[k])
            seg = pboc.phase_segmentation(phase_im, thresh)
            ints = pboc.extract_intensity(seg, fluo_im)

            # Append the image single-cell intensities to our storage vector.
            for val in ints:
                intensities.append(val)

        # Compute the mean for the given strain.
        mean_int = np.mean(intensities)

        # Decide what vector to add this to. It's important to not mix this up.
        if rep_names[j] == 'auto':
            auto_means.append(mean_int)
        if rep_names[j] == 'delta':
            delta_means.append(mean_int)

# Now that we have the two static parts of the numerator and denominator
# determined, let's extract the means and calculate the fold change for the
# RBS mutant strains.

# Iterate through each operator.
for i in range(len(ops)):

    # Set up a storage list for the extracted fold-change values.
    fc = []

    # Iterate through each RBS mutant. Note that we start at '2' here.
    for j in range(2, len(rep_names)):

        # Get the phase and fluorescence image names.
        phase_names = glob.glob('data/lacI_titration/' + ops[i] + '_' +
                                rep_names[j] + '_phase*.tif')
        fluo_names = glob.glob('data/lacI_titration/' + ops[i] + '_' +
                               rep_names[j] + '_yfp*.tif')

        # Set up a storage list for the single-cell intensities.
        intensities = []

        # Perform the image processing.
        for k in range(len(phase_names)):
            phase_im = skimage.io.imread(phase_names[k])
            fluo_im = skimage.io.imread(fluo_names[k])
            seg = pboc.phase_segmentation(phase_im, thresh)
            ints = pboc.extract_intensity(seg, fluo_im)

            # Add the single-cell intensities for a single image to the
            # storage list.
            for val in ints:
                intensities.append(val)

        # Compute the mean intensity for this repressor copy number.
        mean_int = np.mean(intensities)

        # Compute the fold change
        fc.append((mean_int - auto_means[i]) / (delta_means[i] -
                                                auto_means[i]))

    # Add the fold change to the correct storage vector.
    if ops[i] == 'O1':
        O1_fc = fc
    if ops[i] == 'O2':
        O2_fc = fc
    if ops[i] == 'O3':
        O3_fc = fc

# Let's test our theory! We'll plot the experimental data directly ontop of our
# theoretical predictions we made earlier.
plt.plot(rep_numbers, O1_fc, 'ro', label='O1 data')
plt.plot(rep_numbers, O2_fc, 'go', label='O2 data')
plt.plot(rep_numbers, O3_fc, 'bo', label='O3 data')
plt.xlabel('number of repressors')
plt.ylabel('fold-change')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.show()


# Not too bad! While we have a small degree of deviance from the prediction,
# we've still done an impressive job at quantitatitvely describing the behavior
# of the system. Note, however, that the deviance is much larger at the far
# end of the repressor copy numbers (R > 1000). This is because at this high
# of a concentration of repressor, the sensitivity of our microscope actually
# becomes limiting. It's is hard to collect enough photons from YFP when the
# repression is this strong.

# This wraps up our image processing project using python.
