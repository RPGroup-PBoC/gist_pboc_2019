# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Image processing utilities
import skimage.io
import skimage.morphology
import skimage.filters

plt.close('all')
# In this script, we'll learn some basic facts about images, how to read them
# in Python, and the basic principles of thresholding. On the course website,
# (rpgroup-pboc.github.io/gist_pboc_2017), you will find a large image set
# that we will use for the whole project. Before we begin, let's talk a bit
# about the nature of the image set.

# These images are of  E. coli cells with a variety of different copy numbers
# of the LacI repressor molecule. This set is composed of three different LacO
# sequences (O1, O2, and O3), and a variety of different repressor copy numbers
# (indicated by the `R` in the image file name).  In these strains (and the #
# `wt` strain), the LacI repressor molecule represses the expression of a
# Yellow Fluorescent Protein molecule. With more repressor around, less YFP
# molecules are made. There are three strains without an `R` label. These are
# `auto` which is expressing no YFP at all, `delta` which is constitutively
# expressing YFP (has no repressors), and `wt` which has the wild-type number
# of LacI repressors, 22 per cell. Every strain is constitutively expressing
# an mCherry molecule off of a plasmid carried in the cell. This was used to
# make the image segmentation more simple.

# In our project, we will quantify the fold-change in gene expression under
# different repressor copy numbers. To do so, we will need to make measurements
# of the single-cell fluorescence intensities in our images. We'll start by
# learning about images and how to process them in the Python programming
# language.

# It is important to remember that an image is nothing but data -- it is an
# array of points with a specific value. These points are called 'pixels'. The
# values that these pixels can take is related to the construction of the
# camera and is measured as 'bit depth'. To determine the range of pixel
# values in an N bit image ca take, we simply need to compute 2^N - 1. This
# subtraction of 1 is because 0 can be a pixel value as well. For example,
# a 16-bit image can have pixels on the range of 0 -> (2^16 -1 ) = 0 -> 65535.
#  Let's begin by loading an example image into Python.

# Define the image name. We'll look at a constitutively expressing case.
image = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')

# Let's look at the values of the image.
np.shape(image)

# We see that this is a 512 x 512 x 3 image. This means it is 512 pixels wide,
# 512 pixels tall, and has three channels. These three channels correspond to
# the phase contrast image, the mCherry image, and the YFP image. Let's take a
# a look at the image and split it up by the channel.

plt.figure()
plt.imshow(image, cmap=plt.cm.Greys_r)
plt.show()

# We see that bacterial cells are black against a light colored background.
# We can take our 'cross' tool in python and hover over the images to see that
# the pixels with in the bacerium are lower than those of the background. We
# could select only the pixels of the bacterium by drawing a threshold at some
# value and saying anything below that value is 'bacterial'. To figure out
# what this threshold should be, we can look at the image histogram and pick a
# value.

# Generate the image histogram.
plt.figure()
plt.hist(image.flatten(), bins=1000)
plt.xlabel('pixel value')
plt.ylabel('frequency')
plt.show()


# We can clearly see two humps in this image. The left most hump are the pixels
# within our bacterial cells while the major hump are the pixels of the
# background. One could imagine selecting only the bacteria in this image
# by choosing a threshold value between these two peaks and identifying any
# pixel below this threshold as bacterial. Looking at the histogram, we can
# choose a threshold of around 35000 counts.
im_thresh = image < 3500

plt.figure()
plt.imshow(im_thresh, cmap=plt.cm.Greys_r)
plt.show()

# This seems to do a pretty good job at separating the cells from the
# background, but it is really dependent on the actual values recorded by the
# camera. These values can vary from sample to sample (and even within a single
# image). However, the ratio of the cell interior to the background should be
# approximately constant. Let's convert this image to a float (values ranging
# between 0 and 1) and choose this threshold.

im_float = (image - image.min()) / (image.max() - image.min())
plt.figure()
plt.hist(im_float.flatten(), bins=1000)
plt.xlabel('normalized pixel count')
plt.ylabel('counts')
plt.show()

# Our new threshold is somwhere around 0.15. Let's see how that looks.
im_float_thresh = im_float < 0.2
plt.figure()
plt.imshow(im_float_thresh, cmap=plt.cm.Greys_r)
plt.show()


# Why are we getting so much of the background in our segmentation? This is
# because the illumination of the image is not uniform. We can see that the
# lefthand side of the image is darker than that on the right. We can easily
# correct for that by performing a background subtraction. To do so, we will
# very heavily blur the image and subtract it from the original. This will
# remove any large scale aberrations in the intensity of the image leaving
# small scale features (such as bacteria) alone. Let's go a head and give it
# a shot.

# Perform the background subtraction
blur_radius = 50.0  # in units of pixels
im_blur = skimage.filters.gaussian(im_float, sigma=blur_radius)
im_sub = im_float - im_blur

# Now let's look at the image.
plt.figure()
plt.imshow(im_float, cmap=plt.cm.viridis)
plt.title('origiinal')
plt.figure()
plt.show()
plt.imshow(im_sub, cmap=plt.cm.viridis)
plt.title('background subtracted')
plt.show()


# We can see that did the trick. Let's get our new (and final) threshold
# value.
plt.figure()
plt.hist(im_sub.flatten(), bins=1000)
plt.xlabel('normalized pixel value')
plt.ylabel('counts')
plt.show()

# From this histogram, I would choose a threshold value of around -0.2
im_sub_thresh = im_sub < -0.2
plt.figure()
plt.imshow(im_sub_thresh, cmap=plt.cm.Greys_r)
plt.show()

# This seem s very nice, but we are still picking up some garbage in the
# background. We could remove this by only selecting cells that meet a certain
# area threshold, but how can we compute the area of each cell? Right now, all
# pixels that are deemed to be "bacterial" are labeled as 1.0. To the computer
# these are all the same object. We can individually label our segmented cells
# by using the `skimage.measure.label` function, which will individually label
# islands of pixels.
im_lab, num_obj = skimage.measure.label(im_sub_thresh, return_num=True)
print("We've segmented " + str(num_obj) + " objects.")

# Show the labeled image.
plt.figure()
plt.imshow(im_lab, cmap=plt.cm.spectral_r)
plt.colorbar()
plt.show()


# Nice! In the next script, we will learn how to extract some properties of
# each segmented cell and apply an area filter. We can see that even though
# there are only about 24 cells, we've identified 200 objects!
