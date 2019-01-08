# Import our modules

# For doing math
import numpy as np

# For plotting
import matplotlib.pyplot as plt
import seaborn

# For image loading
import skimage.io

# For getting filenames
import glob

a = 2
b = 3 * a
# print(b)

# Load an example image
im = skimage.io.imread('data/ecoli_growth/ecoli_TRITC_18.tif')

# Show our example image
plt.figure()
plt.imshow(im, cmap='gray')
plt.show()

# Generate a histogram of values
im_values = im.flatten()
plt.figure()
plt.hist(im_values, bins=50)
plt.xlabel('pixel value')
plt.ylabel('pixel counts')
plt.show()


# Threshold our image.
im_thresh = im > 207
plt.figure()
plt.imshow(im_thresh, cmap='gray')
plt.show()

# Get the filenames.
filenames = glob.glob('data/ecoli_growth/ecoli_TRITC*.tif')
cell_area = []
for f in filenames:
    im = skimage.io.imread(f)
    im_thresh = im > 207
    area = np.sum(im_thresh)
    cell_area.append(area)

# plot the area as a function of time.
time = np.arange(0, len(filenames), 1) * 5
plt.figure()
plt.plot(time, cell_area, 'o')
plt.xlabel('time (min)')
plt.ylabel('area (pixels)')
plt.show()

# Fit our doubling time.
linear_fit = np.polyfit(time, np.log(cell_area), 1)
slope, intercept = linear_fit

# Calculate the doubling time.
t_double = np.log(2) / slope
print(t_double)


# Plot our fit values.
plt.figure()
fit_curve = intercept + slope * time
plt.plot(time, fit_curve, '-', label='fit')
plt.plot(time, np.log(cell_area), 'o', label='data')
plt.xlabel('time (min)')
plt.ylabel('area (pixels)')
plt.legend()
plt.show()
