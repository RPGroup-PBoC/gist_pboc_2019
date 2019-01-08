# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# For easily grabbing several file names
import glob

# For image processing functions
import skimage.io
import skimage.exposure
import skimage.filters
import skimage.morphology

# In this script, we will compute the growth rate of a bacerial colony over
# time through high-magnification microscopy. The provided images are of
# E. coli cells growing on a hard agar substrate supplemented with a rich
# medium to grow in. The images are at a magnification of 100x and all images
# were taken at 5 minute intervals. In the coming code, we will perform some
# basic segmentation on these images and compute the growth rate.

# We'll start by looking at looking at an image some where in the middle of the
# movie. At this point, there are a bunch of cells stacked next to each other
# forming a microcolony. Our goal will be to figure out some way of measuring
# the amount of area of the image is composed of bacterial cells. To load this
# image, we will pull from an incredibly useful image processing module called
# sci-kit image. We imported this (along with some other modules) above by
# importing skimage.

plt.close('all')
# Load the example image.
im = skimage.io.imread('data/ecoli_growth/ecoli_phase_24.tif')

# Show it
plt.figure()
plt.imshow(im, cmap=plt.cm.Greys_r)
plt.show()

# Noteice that all of teh cells are black and the background is a light gray.
# Remember, an image is data -- it is simply a two-dimensional array of pixel
# values. By hovering over the above image with our mouse, we can see that
# all of the cells are dark (having values of a few hundred counts) while the
# background is much higher. Note that illumination differences between frames
# may be slightly different, making it more difficult to draw a very strict
# boundary for the threshold of what is a cell. However, the approximate ratio
# of the cell pixel values to the background values should be constant. We
# can convert this image to a floating point image which will rescale all of
# the vales from 0 to 1.0

# Make the image a float
im_float = (im - im.min()) / (im.max() - im.min())

# How do we distinguish what is cell mass and what is background? We could
# hover over a bunch of cells to get the idea of what the pixel counts are, but
# a more useful method is to look at the image histogram. This will show us
# the frequency of a given pixel count in the image.

# Show the image histogram.
plt.figure()
plt.hist(im_float.flatten(), bins=500)
plt.xlabel('pixel value')
plt.ylabel('frequency')
plt.show()

# We can see that there are two humps in the histogram. The leftmost hump are
# most likely the cells while the largest hump is likely the  background
# pixels. Lets just show every pixel that is below the value of X to see if
# those are really cells.
plt.figure()
plt.imshow(im_float < 0.27, cmap=plt.cm.Greys_r)
plt.show()


# We did a good job of only selecting the bacteria, however, we also got a
# bunch of the background! This is because the illumination is uneven meaning
# that the left-most part of the image is darker than the right-most part. We
# can correct for this by doing something called a background subtraction. To
# perform this, we'll very heavily blur the image and subtract it from the
# original float image. This means that large variances in illumination will
# be removed while the smaller structures (such as the bacteria) will be
# preserved. Let's blur the image and perform the subtraction.

im_blur = skimage.filters.gaussian(im_float, sigma=30.0)
im_sub = im_float - im_blur
plt.figure()
plt.imshow(im_sub, cmap=plt.cm.Greys_r)
plt.title('background subtracted')
plt.show()


# That is much better! Let's look at the histogram again and try to choose a
# better threshold.
plt.figure()
plt.hist(im_sub.flatten(), bins=500)
plt.xlabel('pixel value')
plt.ylabel('frequency')
plt.title('background subtracted image histogram')


# We see the two peaks much more clearly now. Also, it makes sense that the
# peak of the background pixels is the most frequent as most pixels in our
# image are actually background. Let's go a head and apply the threshold.
im_thresh = im_sub < -0.05

plt.figure()
plt.imshow(im_thresh, cmap=plt.cm.Greys_r)
plt.show()

# We did a much better job, but there are a bunch of small dots around. These
# are stray pixels in the background which are falling below our threshold.
# We can get rid of these by removing all of the objects which are less than
# 50 square pixels. There is a simple command for this packaged in the
# skimage module.
im_large = skimage.morphology.remove_small_objects(im_thresh, min_size=50)
plt.figure()
plt.imshow(im_large, cmap=plt.cm.Greys_r)
plt.show()

# That is much better! But how do we measure area? Remember, by doing the
# thresholding, we are generating a binary image with a pixel value of 1
# wherever there are 'bacterial cells' and 0 else where. In order to get the
# bacterial area of teh image, we can simply sum up all values of the image!
bacterial_area = np.sum(im_large)
print('Bacterial area in this frame is ' + str(bacterial_area) + ' sq pixels')


# In order to measure the growth rate of this colony, we will want to
# iterate across all images in this experiment, repeat the steps above, and
# calculate the area. We can then plot the bacterial area as a function of
# time and get a measure of the growth rate.

# We'll start by getting a list of all of the image names in the file. We can
# do this using another python module called glob.
image_names = glob.glob('data/ecoli_growth/ecoli_phase_*.tif')

# The asterisk means it will get all file names that match that pattern where
# anything can occur betweeen `ecoli_phase_ ` and `.tif`. With this set of
# file names, now we can simply iterate through each file and perform the same
# set of steps!
cell_area = np.zeros(len(image_names))  # Make an empty storage vector.
for i in range(len(image_names)):
    # Load the image.
    im = skimage.io.imread(image_names[i])

    # Make the image a float.
    im_float = (im - im.min()) / (im.max() - im.min())

    # Perform the background subtraction.
    im_blur = skimage.filters.gaussian(im_float, sigma=30.0)
    im_sub = im_float - im_blur

    # Apply the threshold
    im_thresh = im_sub < -0.05
    im_large = skimage.morphology.remove_small_objects(im_thresh, min_size=10)

    # Compute the cell area.
    cell_area[i] = np.sum(im_large)


# Let's plot the bacterial area as a function of time.
time = np.arange(0, len(image_names) * 5, 5)
plt.figure()
plt.plot(time, cell_area, 'o')
plt.xlabel('time (min)')
plt.ylabel('cell area (sq. pixels)')
plt.show()

# Since we predict that the growth is exponential with time, this plot should
# be approximately linear on a log-Y scale.
plt.figure()
plt.plot(time, np.log(cell_area), 'o')
plt.xlabel('time (min)')
plt.ylabel('log(cell_area (sq. pixels))')
plt.show()

# That is impressively linear given how rough our segmentation algorithm is.
# Let's fit this trend to an exponential curve to find the doubling time. If
# we assume that this growth is exponential, then we can fit
#
#       A_t = A_0 * exp(k * t)
#       ln(A_t) = ln(A_0) + k * t.
#
# To determine the doubling time, this is a simple rearrangement of the above
# linear equation to
#
#       t_double = ln(2) / k.
#
# To do this, we'll use the NumPy polyfit function, although there are many
# different available utilities to do this.
linear_fit = np.polyfit(time, np.log(cell_area), 1)

# The output of this function is an array with the  slope an intercept.
slope, intercept = linear_fit

# Now, let's compute the doubling time.
t_double = np.log(2) / slope
print('The doubling time is ' + str(t_double) + ' min.')


# Let's also plot the fit on our raw data.
plt.figure()
fit_curve = intercept + slope * time
plt.plot(time, fit_curve, 'k-', label='fit')
plt.plot(time, np.log(cell_area), 'o', label='experiment')
plt.xlabel('time (min)')
plt.ylabel('log(cell area (sq. pixels))')
plt.show()


# That is a pretty good fit! The compute value for the doubling time also makes
# some sense with our intuition of how cells growth. In class, we determined
# that our rule of thumb for E. coli growth is around twenty minutes whereas
# our fitted value is a little bit over that. There are a few reasons for this.
# First, our segmentation isn't 'perfect'. We aren't really segmenting most of
# the cells when the colony is fully developed. Secondly, these cells are
# growing in an aerobic environment sandwiched between the pad  and the glass
# under the microscope. However, in just around 60 lines of acutal code, we
# are able to get a pretty good measure!
