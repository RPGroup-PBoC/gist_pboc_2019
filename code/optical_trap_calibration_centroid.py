# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# For image processing
import skimage.io
import skimage.filters
import skimage.segmentation
import skimage.measure

# In this script, we will calibrate the strength of an optical trap. The data
# set for this script is of a one micron bead trapped in a laser optical trap
# with a 5.2x beam expansion from a 1 mm diameter red laser source. This image
# was taken at 30 frames per second and the interpixel distance of the camera
# is 42 nm per pixel at 22 C. As we discussed in class, we know that the force of
# the trap is related to the mean squared displacement of the bead by
#
#   F_trap = <(x - x_avg)^2> / kBT
#
# where <(x - x_avg)^2> is the mean squared displacement of the bead, kB is the
# Boltzmann constant, and T is the temperature of the system. Our goal is to
# determine the mean squared displacement of the bead in the trap. To do so,
# we will use some image processing techniques to segment the bead and identify
# the centroid. You should note that if we were calibrating this optical trap
# for "real-life" experiments, we would use more sophisticated techniques to
# calculate the trap to a sufficient precision.

# To begin, let's load up the time series of our bead and look at the first
# image.
bead_ims = skimage.io.imread('data/optical_tweezer/trapped_bead_5.2x_4_MMStack_Pos0.ome.tif')
plt.figure()
plt.imshow(bead_ims[0], cmap=plt.cm.Greys_r)
plt.show()

# We see that the bead is dark on a light background with some fringing on
# the side. We took these images such that the bead was dark to simplify our
# segmentation analysis.

# Because these images are relatively clean, we can imagine using a threshold
# to determine what is bead and what is not. However, there is a lot of noise
# in the background of this image. Before we identify a threshold, let's try
# blurring the image with a gaussian blur to smooth it out. We'll then look
# at the histogram.

im_blur = skimage.filters.gaussian(bead_ims[0], sigma=1)
plt.figure()
plt.imshow(im_blur, cmap=plt.cm.Greys_r)
plt.show()

# Now, let's look at the image histogram to choose a threshold.
plt.figure()
plt.hist(im_blur.flatten(), bins=1000)
plt.xlabel('pixel value')
plt.ylabel('count')
plt.show()

# It seems pretty distinct what pixels correlate to our bead. Let's impose a
# threshold value of 0.2 and see how well it works.
threshold = 0.15
im_thresh = im_blur < threshold
plt.figure()
plt.imshow(im_thresh, cmap=plt.cm.Greys_r)
plt.show()

# That seems to do a pretty good job! However, there is an object that is touching the border of the image. Since we only want to get the bead segmented, we'll clear the segmentation of anything that is touching the border.
im_border = skimage.segmentation.clear_border(im_thresh)

# Now, we want to find the position of the middle of the bead. While it would
# be best to find the center to sub-pixel accuracy by fitting a two-dimensional
# gaussian, we'll find it by extracting the centroid from the regionprops
# measurement. We'll first label the image and then extract the properties.
im_label = skimage.measure.label(im_border)
props = skimage.measure.regionprops(im_label)
centroid_x, centroid_y = props[0].centroid

# Now, let's plot the position of the centroid on the image of our bead to see
# how well it worked. Remember that images are plotted y vs x, so we will need
# to plot the centroid_y first.
plt.figure()
plt.imshow(im_blur)
plt.plot(centroid_y, centroid_x, 'o')
plt.show()

# That's pretty good! Now, to determine the mean squared  displacement, we
# will want to find the position of the bead at each frame. To do this, we'll
# loop through all of the images and perform the exact same operations.
centroid_x, centroid_y = [], []  # Empty storage lists.
for i in range(len(bead_ims)):
    # Blur the image.
    im_blur = skimage.filters.gaussian(bead_ims[i], sigma=1)

    # Segment the image.
    im_thresh = im_blur < threshold

    # Clear the border.
    im_border = skimage.segmentation.clear_border(im_thresh)
    im_large = skimage.morphology.remove_small_objects(im_border)

    # Label the image and extract the centroid.
    im_lab = skimage.measure.label(im_large)
    props = skimage.measure.regionprops(im_lab, intensity_image=np.invert(bead_ims[i]))
    x, y = props[0].weighted_centroid

    # Store the x and y centroid positions in the storage list.
    centroid_x.append(x)
    centroid_y.append(y)

# Now, let's generate some plots to see if our analysis makes sense. We'll
# plot the centroid x vs centroid y position to make sure that our bead seems
# to be diffusing in the expected manner. We'll also plot the x and y positions
# as a function to time to see if there are any bumps in the table or if our
# segmentation went awry.
plt.figure()
plt.plot(centroid_x, centroid_y, '-')
plt.xlabel('x position (pixels)')
plt.ylabel('y position (pixels)')

# Now plot them as a function of time.
time_vec = np.arange(0, len(bead_ims), 1) #* (1 / 50)  # Converted to seconds.
plt.figure()
plt.plot(time_vec, centroid_x, '-')
plt.xlabel('time (s)')
plt.ylabel('x position (pixels)')

plt.figure()
plt.plot(time_vec, centroid_y, '-')
plt.xlabel('time (s)')
plt.ylabel('y position (pixels)')
plt.show()
# It looks like something gets bumped at frame 24 and then around 120. Let's
# restrict our analysis to that zone. That all looks good! It seems to be
# diffusing as expected. Now let's # calculate the mean squared displacement
# and compute the trap force.
ip_dist = 0.042  # Physical distance in units of  microns per pixel
centroid_x_micron = np.array(centroid_x) * ip_dist
centroid_y_micron = np.array(centroid_y) * ip_dist

# Compute the means and msd.
mean_x = np.mean(centroid_x_micron[24:120])
mean_y = np.mean(centroid_y_micron[24:120])
msd_x = np.mean((centroid_x_micron[24:120] - mean_x)**2)
msd_y = np.mean((centroid_y_micron[24:120] - mean_y)**2)

# Compute the trap force.
kT = 4.1E-3  # In units of pN * micron
k_x = kT / msd_x
k_y = kT / msd_y
print('Trap force in x dimension is ' + str(k_x) + ' pN micron')
print('Trap force in y dimension is ' + str(k_y) + ' pN micron')

# Wow! That's a strong trap. Note that this is an approximation of the trap
# force. To precisely measure it, we should determine the centroid of the bead
# to sub-pixel accuracy.
