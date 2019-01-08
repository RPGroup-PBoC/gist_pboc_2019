# Import the necessary modules.
import numpy as np
import skimage.io
import glob
import scipy.ndimage
import skimage.morphology
import skimage.segmentation
import matplotlib.pyplot as plt

# Load the images.
im_glob = glob.glob('data/optical_tweezer/bead*.tif')
x = 'data/optical_tweezer/trapped_bead_5.2x_4_MMStack_Pos0.ome.tif'
im = skimage.io.ImageCollection(im_glob)


disk_radius = 3 # though choose what you like
selem = skimage.morphology.disk(disk_radius)

def center_of_mass(im, selem):

    # Get peak i and j (if multiple maxima, just take first)
    i, j = np.where(im == im.max())
    i = i[0]
    j = j[0]

    # Get the i and j extent (radii) of the structuring element
    r = np.array(selem.shape) // 2
    r_i, r_j = r

    # Get indices of non-zero entries in structuring element for convenience
    ii, jj = np.nonzero(selem)

    # Define indices such that index zero is in center of selem
    i_pos = ii - r_i
    j_pos = jj - r_j

    # Width of structuring element
    w = 2 * r + 1

    # Make subimage that has selem
    sub_im = im[i - r_i:i + r_i + 1, j - r_j:j + r_j + 1]

    # Compute center of mass
    eps_i, eps_j = np.array([np.dot(i_pos, sub_im[ii,jj]), np.dot(j_pos, sub_im[ii,jj])]) / sub_im[ii,jj].sum()

    # Return center of mass of bead
    return i + eps_i, j + eps_j

# find the center of mass
x_pos, y_pos = [], []
for i in range(100):
    im_blur = skimage.filters.gaussian(im[i], 3)
    m = skimage.measure.moments(im[i])
    x = m[0, 1] / m[0, 0]
    y = m[1, 0] / m[0, 0]
    x_pos.append(x)
    y_pos.append(y)


plt.figure()
plt.plot(x_pos, y_pos, 'o')
plt.show()


ip_dist = 0.042  # Physical distance in units of microns per pixel
centroid_x_micron = np.array(x_pos) * ip_dist
centroid_y_micron = np.array(y_pos) * ip_dist

# Compute the means and msd.
mean_x = np.mean(centroid_x_micron)
mean_y = np.mean(centroid_y_micron)
msd_x = np.mean((centroid_x_micron - mean_x)**2)
msd_y = np.mean((centroid_y_micron - mean_y)**2)

# Compute the trap force.
kT = 4.1E-3  # In units of pN * micron
k_x = kT / msd_x
k_y = kT / msd_y
print('Trap force in x dimension is ' + str(k_x) + ' pN micron')
print('Trap force in y dimension is ' + str(k_y) + ' pN micron')
