# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# Image processing utilities.
import skimage.io

# For peak detection
import scipy.signal

plt.close('all')
# This script serves as an example of how one would more accurately
# determine the interpixel distance of a camera. This was taken on an optical
# tweezer set up built by Nathan Belliveau and was captured on a qimaging
# camera. This is a chrome-plated graticule where each division is 10 microns
# apart.

# Load the image.
grat_im = skimage.io.imread('data/optical_tweezer_old/graticule.tif')

# Plot all of the line profiles.
im_shape = np.shape(grat_im)
x_vec = np.arange(0, im_shape[1], 1)
plt.figure()
for i in range(im_shape[0]):
    plt.plot(x_vec, grat_im[0,:], 'k-', alpha=0.5)

# Plot the mean
mean_prof = np.mean(grat_im, axis=0)
plt.plot(x_vec, mean_prof, 'r-')
plt.xlabel('x position')
plt.ylabel('pixel value')
plt.show()

# Find the local min.
local_min = scipy.signal.argrelmin(mean_prof, order=100)[0]

# Plot the points of the minimum for good measure.
for i in range(len(local_min)):
    plt.plot(local_min[i], mean_prof[local_min[i]], 'o')

plt.show()


# Find the distance between each one.
diff_dist = 10.0 / np.diff(local_min)
ip_dist = np.mean(diff_dist)
print(ip_dist)

# This yields a distance of about 42 nm per pixel
