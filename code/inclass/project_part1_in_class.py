# Import our favorite modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load image processing utilities
import skimage.io
import skimage.filters
import skimage.morphology
import skimage.segmentation

# Load an example image
im = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')

# Show the image
plt.figure()
plt.imshow(im, cmap='gray')
plt.show()

# Generate a histogram of the image
im_values = im.flatten()
plt.figure()
plt.hist(im_values, bins=1000)
plt.xlabel('pixel value')
plt.ylabel('counts')
plt.show()

# Threshold the image.
plt.close('all')
im_thresh = im < 2500
plt.figure()
plt.imshow(im_thresh, cmap='gray')
plt.show()

# Normalize the image and plot the histogram.
im_norm = (im - im.min()) / (im.max() - im.min())

im_norm_values = im_norm.flatten()
plt.figure()
plt.hist(im_norm_values, bins=1000)
plt.xlabel('normalized pixel values')
plt.ylabel('frequency')
plt.show()

# Thresh our normalized image
im_norm_thresh = im_norm < 0.3
plt.figure()
plt.imshow(im_norm_thresh, cmap='gray')
plt.show()

plt.close('all')

# Do a background subtraction.
sigma = 50.0
im_blur = skimage.filters.gaussian(im_norm, sigma)

# Show the blur
plt.figure()
plt.imshow(im_blur, cmap='gray')
plt.show()

# Do the subtraction
im_sub = im_norm - im_blur

# Plot it
plt.figure()
plt.imshow(im_sub, cmap='gray')
plt.show()

# Generate the histogram of the flattened image
im_sub_values = im_sub.flatten()
plt.figure()
plt.hist(im_sub_values, bins=1000)
plt.xlabel('flattened pixel value')
plt.ylabel('counts')
plt.show()

# Apply a new threshold value
im_sub_thresh = im_sub < -0.2
plt.figure()
plt.imshow(im_sub_thresh, cmap='gray')
plt.show()

# Remove the small background pixels.
im_large = skimage.morphology.remove_small_objects(im_sub_thresh, min_size=50)
plt.figure()
plt.imshow(im_large, cmap='gray')
plt.show()
