# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Image processing utilities
import skimage.io
import skimage.filters
import skimage.segmentation
import skimage.measure

# In this script, we will learn some more principles regarding image processing
# and segmentation. We'll examine how to extract properties from segmented
# objects and show off our fancy segmentation mask by making an overlay.

# Before we begin, let's type what we did in the first step of the project.
image = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')

# Do our normalization, segmentation, and labeling.
im_float = (image - image.min()) / (image.max() - image.min())
im_blur = skimage.filters.gaussian(im_float, sigma=50.0)
im_sub = im_float - im_blur
im_thresh = im_sub < -0.2
im_lab, num_cells = skimage.measure.label(im_thresh, return_num=True)

plt.figure()
plt.imshow(im_lab, cmap=plt.cm.spectral)
plt.show()

# Looking at the above image, we see that there are far more segemented objects
# than there are actual cells. This is because we are segmenting some of the
# pixels in the background of the images. We imagined that we could get rid of
# these pixels by slecting objects wihch meet a set of area bounds. Before we
# apply any bounds, let's just look at the areas of all of the cells in our
# image.

# Set the physical distance of the pixels in our camera
ip_dist = 0.160  # in units of um/pixel

# Set up a list where we will store the values of the comptued cell areas.
area = []

# We can get the properties of each region in our image by using
# skimage.measure.regionprops. This returns properties such as area, label,
# mean intesity, eccentricity, image moments,  etc.
props = skimage.measure.regionprops(im_lab)

# We'll now iterate through each property and extract the area.
for prop in props:
    area.append(prop.area * ip_dist**2)


# Let's take a look at the distribution of cell areas.
plt.figure()
plt.hist(area, bins=75)
plt.xlabel('object area (sq. micron)')
plt.ylabel('counts')
plt.show()


# Yikes, it seems like we have a bunch of garbage. What would some good bounds
# be? Our rule-of-thumb is that E. coli is about 2 microns long by one micron
# wide. If we approximate our cell as a rectangle, this gets us to an area of
# 2 sq micron. Of course, not all of our cells are ideal. We can be a little
# more lenient and say that our smallest cell would probably be 0.5 sq micron
# with our largest being about 6 sq. micron. We see in our histogram that we
# have some distribution between  about 1.5 - 3.5 square micron. There is
# another distribution that is much smaller than our bounds. Let's filter out
# those objects and see what we segment.

approved_obj = np.zeros_like(im_lab)
for prop in props:
    obj_area = prop.area * ip_dist**2
    if (obj_area > 0.5) & (obj_area < 6):
        # This is an approved cell, let's store the label.
        approved_obj += (im_lab == prop.label)


plt.figure()
plt.imshow(approved_obj, cmap=plt.cm.spectral)
plt.show()

# That looks pretty good! Let's make sure this makes sense by looking at the
# new histogram. This means we will have to relabel the image.
im_relab, num_cells = skimage.measure.label(approved_obj, return_num=True)
print("We've identified " + str(num_cells) + " cells!")

cell_props = skimage.measure.regionprops(im_relab)
cell_areas = []
for prop in cell_props:
    cell_areas.append(prop.area * ip_dist**2)

plt.figure()
plt.hist(cell_areas, bins=10)
plt.xlabel('cell area (sq. micron)')
plt.ylabel('counts')
plt.show()

# That looks great! We've even selected the right number of cells. We have one
# last issue, though. There seems to be a cell on the edge that is not
# completely in the image. Since we are ultimately interested in extracting
# quantitative information about the fluorescence, we want to remove any cells
# that are not completely in the field of view. To do this, we can use the
# skimage.segmentation.clear_border command which will delete all cells that
# are touching the border. This only works on a binary image, meaning that we
# will have to relabel the image after the cell has been removed.

im_border = skimage.segmentation.clear_border(approved_obj)
im_border_lab = skimage.measure.label(im_border)
plt.figure()
plt.imshow(im_border_lab, cmap=plt.cm.spectral_r)
plt.show()

# Very nice! We now have a segmentation mask that passes our tests. We've done
# all of our segmentation in phase contrast, which is not the chanel which has
# our fluorescence information. Let's load up that image and take a look.
im_fluo = skimage.io.imread('data/lacI_titration/O2_delta_yfp_pos_16.tif')
plt.figure()
plt.imshow(im_fluo, cmap=plt.cm.Greys_r)
plt.show()

# This is a much different image. The cells are now bright against dark
# background. Since this is the channel we are interested in, why didn't we
# segment this image? Well, we are interested in the quantitative information.
# If we segmented through simple thresholding (like we've been doing), we would
# preferentially segment the brightest cells while ignoring the cells that may
# have very little to no signal. We can use the mask we generated from
# segmenting in phase contrast and apply it on this image to extract the
# fluorescence information.

# Extract the fluorescence properties.
cell_fl_props = skimage.measure.regionprops(im_border_lab,
                                            intensity_image=im_fluo)

# Let's look at the distribution of mean cell intensities.
cell_ints = []
for prop in cell_fl_props:
    cell_ints.append(prop.mean_intensity)

plt.figure()
plt.hist(cell_ints, bins=10)
plt.xlabel('fluorescence pixel intensity (a. u.)')
plt.ylabel('counts')
plt.show()

# We can see that the fluorescence intensity is pretty bright and distributed
# between 2000 and 6000 counts for 23 cells in our segmentation mask. We can
# very easily compute the mean intensity of the cells in this image using
# the list of cell intensities we just generated.
mean_intensity = np.mean(cell_ints)
print('The mean cell intensity is  ' + str(mean_intensity) + ' counts.')


# Let's do one more thing for fun. It's often useful to generate an overlay
# of your segmentation mask to show others that you are only looking at cell
# mass in your measurements. We'll do this by generating an image where the
# segmented cells are colored in blue over our original image.

# Make a copy of our float phase image.
phase_copy = np.copy(im_float)
# Color the segemented parts.
phase_copy[im_border > 0] = 1.0

# Make an RGB image with the blue channel phase_copy.
merge = np.dstack((im_float, im_float, phase_copy))

# Show it!
plt.figure()
plt.imshow(merge)
plt.show()

# Looks great! We've done a lot of work in this tutorial. In the next script,
# you will get the chance to write everything we've done here as functions so # it is easy to apply it to a large stack of images.
