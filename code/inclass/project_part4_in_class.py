# By importing stuff.
import numpy as np
import matplotlib.pyplot as plt
import seaborn
seaborn.set_context('talk')
# To load images.
import glob
import skimage.io

# Course utilities
import pboc_utils as pboc

# Test our theory.
R = np.logspace(0, 4)  # Number of repressors per cell
Nns = 5E6  # Number of nonspecific binding sites
ep_r = -13.9  # in units of kT

# Compute the fold-change
fold_change_prediction = 1 / (1 + (R / Nns) * np.exp(-ep_r))
# Plot our prediction.
plt.figure()
plt.plot(R, fold_change_prediction, label='theory')
plt.xlabel('number of repressors')
plt.ylabel('fold-change')
plt.legend()
plt.ylim([0.001, 1.1])
plt.xscale('log')
plt.yscale('log')
plt.show()

#  Get the auto and delta sets first.
auto_phase_files = glob.glob('../data/lacI_titration/O2_auto_phase*.tif')
auto_yfp_files = glob.glob('../data/lacI_titration/O2_auto_yfp*.tif')
auto_intensities = []
thresh = -0.2
for i in range(len(auto_phase_files)):
    # Load the phase image.
    im = skimage.io.imread(auto_phase_files[i])
    yfp_im = skimage.io.imread(auto_yfp_files[i])

    # Do the segmentation
    mask = pboc.phase_segmentation(im, thresh)
    ints = pboc.extract_intensity(mask, yfp_im)
    for value in ints:
        auto_intensities.append(value)

# Compute the mean auto intensity.
mean_auto = np.mean(auto_intensities)

# Load the delta files.
delta_phase_files = glob.glob('../data/lacI_titration/O2_delta_phase*.tif')
delta_yfp_files = glob.glob('../data/lacI_titration/O2_delta_yfp*.tif')
delta_intensities = []
for i in range(len(delta_phase_files)):
    # Load the phase image.
    im = skimage.io.imread(delta_phase_files[i])
    yfp_im = skimage.io.imread(delta_yfp_files[i])

    # Do the segmentation
    mask = pboc.phase_segmentation(im, thresh)
    ints = pboc.extract_intensity(mask, yfp_im)
    for value in ints:
        delta_intensities.append(value)

# Compute the mean delta intensity.
mean_delta = np.mean(delta_intensities)


# Loop over all O2 images.
rep_names = ['wt', 'R60', 'R124', 'R260', 'R1220', 'R1740']
rep_number = [11, 60, 124, 260, 1220, 1740]
fold_change = []
for i in range(len(rep_names)):
    # Load the file names
    phase_files = glob.glob('../data/lacI_titration/O2_' + rep_names[i] + '_phase*.tif')
    yfp_files = glob.glob('../data/lacI_titration/O2_' + rep_names[i] + '_yfp*.tif')

    # Extract the mean intensities
    sample_intensities = []
    for j in range(len(phase_files)):
        # Load the phase image.
        im = skimage.io.imread(phase_files[j])
        yfp_im = skimage.io.imread(yfp_files[j])

        # Do the segmentation
        mask = pboc.phase_segmentation(im, thresh)
        ints = pboc.extract_intensity(mask, yfp_im)
        for value in ints:
            sample_intensities.append(value)

    # Compute the mean delta intensity.
    mean_int = np.mean(sample_intensities)

    # Calculate fold change.
    fold_change_experiment = (mean_int - mean_auto) / (mean_delta - mean_auto)

    # Save the value.
    fold_change.append(fold_change_experiment)

# Plot our data!!
plt.plot(rep_number, fold_change, 'o', label='experiment')
plt.legend()
plt.show()
plt.tight_layout()
plt.savefig('../figures/theory_experiment_agreement.png', bbox_inches='tight')
