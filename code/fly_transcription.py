# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# For getting file names.
import glob

# For image processsing.
import skimage.io
import skimage.filters
import skimage.segmentation
import skimage.measure

# In this exercise, we will measure the rate of transcription in a nuclear
# cycle 14 of a developing Drosophila melanogaster embryo. Specificially, we
# are watching the beginning and end of transcription of the Hunchback
# morphogen. This was done in a particularly clever way by using the MS2 system
# to tag the 5' and 3' ends of the mRNA in separate embryos. In this
# experimental systems, fluorescent puncta appear when that region of the mRNA
# is properly transcribed. By watching the appearance of the relevant spots, we
# can make a measurement of the rate of transcription and investigate whether
# the timing makes sense, given the time in between nuclear divisions. It is
# important to know that the length of this gene is 3361 base pairs.

# In this exercise, we will measure the rate of transcription in a nuclear
# cycle 14 of a developing Drosophila melanogaster embryo. Specificially, we
# are watching the beginning and end of transcription of the Hunchback
# morphogen. This was done in a particularly clever way by using the MS2 system
# to tag the 5' and 3' ends of the mRNA in separate embryos. In this
# experimental systems, fluorescent puncta appear when that region of the mRNA
# is properly transcribed. By watching the appearance of the relevant spots, we
# can make a measurement of the rate of transcription and investigate whether
# the timing makes sense, given the time in between nuclear divisions.

# To begin with, let's take a look at an image somewhere in the middle of the
# movie for the 5' labeled mRNA case.
im = skimage.io.imread('data/fly_elongation/5Loops0200.tif')
plt.figure()
plt.imshow(im, cmap=plt.cm.viridis)
plt.show()

# How can we count? It seems like we would be able to threshold
# the images to find the spots. To pick a threshold, let's take a look at the
# histogram.
plt.figure()
plt.hist(im.flatten(), bins=1000)
plt.xlabel('pixel value')
plt.ylabel('counts')

# We see that we have a very large peak which corresponds to the background
# pixels and a long tail to the right which likely corresponds to our spots.
# Let's try applying a threshold of 0.11 and see how we've done.
im_thresh = im > 900
plt.figure()
plt.imshow(im_thresh, cmap=plt.cm.viridis)
plt.show()

# That seems to do a pretty good job! But we still have some spots right on the
# edge which are unclear if they are spots or not. To prevent our selves from
# overcounting, we can remove any spots touching the border.
im_border = skimage.segmentation.clear_border(im_thresh)

# We can individually label each spot and count them to double check our
# segmentation.
im_lab, num_spots = skimage.measure.label(im_border, return_num=True)
print('There are ' + str(num_spots) + ' mRNA spots in this image.')


# That is pretty good, it matches our count by eye exactly!  Our goal is to
# count the number of spots per frame and plot their arrival as a function of
# time. With about 300 frames, this is too much to do by hand. Let's write a
# function so this can be automated.
def spot_counter(im, threshold):
    """Counts the number of puncta in an image"""

    # Threshold the image.
    im_thresh = im > threshold

    # Clear the border.
    im_border = skimage.segmentation.clear_border(im_thresh)

    # Label the image and return the number of spots.
    im_lab, num_spots = skimage.measure.label(im_border, return_num=True)
    return num_spots


# Now let's iterate over all of the 5' labeled images. We can get all of the
# files through using glob.
spots_5p_files = glob.glob('data/fly_elongation/5Loops*.tif')
spots_5p = []
for i in range(len(spots_5p_files)):
    im = skimage.io.imread(spots_5p_files[i])
    num_spots = spot_counter(im, 900)
    spots_5p.append(num_spots)

# While we are at it, we can do the same procedure with the 3' labeled mRNA.
spots_3p_files = glob.glob('data/fly_elongation/3Loops*.tif')
spots_3p = []
for i in range(len(spots_3p_files)):
    im = skimage.io.imread(spots_3p_files[i])
    num_spots = spot_counter(im, 900)
    spots_3p.append(num_spots)

# Let's plot the number of spots as a function of time.
time_5p = np.arange(0, len(spots_5p), 1) * 10  # conversion to seconds
time_3p = np.arange(0, len(spots_3p), 1) * 10  # conversion to seconds
plt.figure()
plt.plot(time_5p, spots_5p, '-', label="5' labeled")
plt.plot(time_3p, spots_3p, '-', label="3' labeled")
plt.xlabel('time (s)')
plt.ylabel('spot number')
plt.legend()
plt.show()

# This looks good, but we know that these two image sets are not synchronized.
# By this, I mean that each image does not correspond to the same point in
# development between the two labeling types. By looking at the image stacks
# in ImageJ, we can see that the anaphase of the last nuclear cycle (nuclear
# cycle 14) for the 5' labeled sample begins at frame 99 and at frame 125
# for the 3' labeled sample. We can plot only these ranges to determine the
# rate of transcription.
spots_5p_sync = spots_5p[98:]
time_5p_sync = time_5p[98:]
spots_3p_sync = spots_3p[124:]
time_3p_sync = time_3p[124:]

# Now let's take a look.
plt.figure()
plt.plot(time_5p_sync, spots_5p_sync, label="5' labeled")
plt.plot(time_3p_sync, spots_3p_sync, label="3' labeled")
plt.xlabel('time (s)')
plt.ylabel('spot number')
plt.legend(loc='upper left')
plt.show()

# We can see that the difference in the onset of peaks is about 250 - 300
# seconds. Knowing that the length of the gene is 3361 base pairs, we
# can get to a transcription rate of about 10 - 20 base pairs per second.
# Since a nuclei takes about ten minutes to duplicate, half of that time is
# being spent just transcribing this gene! We haven't even mentioned the
# translation!.
