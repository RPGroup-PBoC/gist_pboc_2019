# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.optimize
import glob
import skimage.io
import skimage.morphology
import scipy.constants
# Define functions from Justin Bois
# Fit symmetric Gaussian to x, y, z data
def fit_gaussian(x, y, z):
    """
    Fits symmetric Gaussian to x, y, z.

    Fit func: z = a * exp(-((x - x_0)**2 + (y - y_0)**2) / (2 * sigma**2))

    Returns: p = [a, x_0, y_0, sigma]
    """

    def sym_gaussian(p):
        """
        Returns a Gaussian function:
        a**2 * exp(-((x - x_0)**2 + (y - y_0)**2) / (2 * sigma**2))
        p = [a, x_0, y_0, sigma]
        """
        a, x_0, y_0, sigma = p
        return a**2 \
                * np.exp(-((x - x_0)**2 + (y - y_0)**2) / (2.0 * sigma**2))

    def sym_gaussian_resids(p):
        """Residuals to be sent into leastsq"""
        return z - sym_gaussian(p)

    def guess_fit_gaussian():
        """
        return a, x_0, y_0, and sigma based on computing moments of data
        """
        a = z.max()

        # Compute moments
        total = z.sum()
        x_0 = np.dot(x, z) / total
        y_0 = np.dot(y, z) / total

        # Approximate sigmas
        sigma_x = np.dot(x**2, z) / total
        sigma_y = np.dot(y**2, z) / total
        sigma = np.sqrt(sigma_x * sigma_y)

        # Return guess
        return (a, x_0, y_0, sigma)

    # Get guess
    p0 = guess_fit_gaussian()

    # Perform optimization using nonlinear least squares
    popt, junk_output, info_dict, mesg, ier = \
            scipy.optimize.leastsq(sym_gaussian_resids, p0, full_output=True)

    # Check to make sure leastsq was successful.  If not, return centroid
    # estimate.
    if ier in (1, 2, 3, 4):
        return (popt[0]**2, popt[1], popt[2], popt[3])
    else:
        return p0

def bead_position_pix(im, selem):
    """
    Determines the position of bead in image in units of pixels with
    subpixel accuracy.
    """
    # The x, y coordinates of pixels are nonzero values in selem
    y, x = np.nonzero(selem)
    x = x - selem.shape[1] // 2
    y = y - selem.shape[0] // 2

    # Find the center of the bead to pixel accuracy
    peak_flat_ind = np.argmax(im)
    peak_j = peak_flat_ind % im.shape[0]
    peak_i = (peak_flat_ind - peak_j) // im.shape[1]

    # Define local neighborhood
    irange = (peak_i - selem.shape[0] // 2, peak_i + selem.shape[0] // 2 + 1)
    jrange = (peak_j - selem.shape[1] // 2, peak_j + selem.shape[1] // 2 + 1)

    # Get values of the image in local neighborhood
    z = im[irange[0]:irange[1], jrange[0]:jrange[1]][selem.astype(np.bool)]

    # Fit Gaussian
    a, j_subpix, i_subpix, sigma = fit_gaussian(x, y, z)

    # Return x-y position
    return np.array([peak_i + i_subpix, peak_j + j_subpix])


# Load the images.
g = 'data/optical_tweezer/trapped_bead_5.2x_4_MMStack_Pos0.ome.tif'
im = skimage.io.imread(g)


# We will use the nine-point estimate (as is typically done)
selem = skimage.morphology.square(3)

# Loop through and find centers
centers = []
length=100
time = np.arange(0, length, 1)
for i in range(length):
    centers.append(bead_position_pix(np.invert(im[i]), selem))

# Store as NumPy array
centers = np.array(centers)

# Get displacements
x = centers[:,1] - centers[:,1].mean()
y = centers[:,0] - centers[:,0].mean()

# Plot displacement
plt.figure()
plt.plot(time, centers[:,0], lw=1, zorder=1, label=r'$x$')
plt.figure()
plt.plot(time, centers[:,1], lw=0.5, zorder=0, label=r'$y$')
plt.xlabel('time (s)')
plt.ylabel('$x$, $y$ (pixels)')
plt.legend(loc='lower left')


# Get x and y in real units
ip_dist = 0.042
x_micron = x * ip_dist
y_micron = y * ip_dist

# Get k's from equipartition
kT = scipy.constants.k * (273.15 + 22.0) * 1e18
k_x = kT / (x_micron**2).mean()
k_y = kT / (y_micron**2).mean()

# Print result
print('k_x = %.2f pN/µm' % k_x)
print('k_y = %.2f pN/µm' % k_y)

plt.figure()
plt.plot(centers[:,0], centers[:,1], '-')
plt.show()
