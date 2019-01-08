# Import our favorite modules
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# For image processing...
import skimage.io
import skimage.measure

# For getting the file names
import glob

# Look at an example image.
im = skimage.io.imread('data/garcia_et_al_2013/3Loops0200.tif')

# Show our example image.
plt.close('all')
plt.figure()
plt.imshow(im, cmap='gray')
plt.show()

# Show a histogram of the image
im_values = im.flatten()
plt.figure()
plt.hist(im_values, bins=700)
plt.xlabel('pixel value')
plt.ylabel('counts')
plt.show()

# Apply a threshold on our image of spots.
im_thresh = im > 1000
plt.figure()
plt.imshow(im_thresh, cmap='gray')
plt.show()

# Label and count the spots.
im_label, num_spots = skimage.measure.label(im_thresh, return_num=True)
print('We found ' + str(num_spots) + ' spots! Does this make sense?')


# Write a function to count spots.
def spot_counter(image, threshold):
    """
    This function counts the number of spots in a given image.
    The image is thresholded with the spots being brighter
    than the threshold value.
    """
    # Do the steps of our segmentation.
    im_thresh = image > threshold

    # Label our image and find the number of spots.
    im_label, num_spots = skimage.measure.label(im_thresh, return_num=True)

    # Now return the value.
    return num_spots

# Load our image names
images_5p = glob.glob('data/garcia_et_al_2013/5Loops*.tif')
spots_5p = []
threshold = 1000
# Loop through and find the spot number.
for f in images_5p:
    # Load the image.
    im = skimage.io.imread(f)

    # Feed it to our function.
    num_spots = spot_counter(im, threshold)

    #Store num_spots
    spots_5p.append(num_spots)

images_3p = glob.glob('data/garcia_et_al_2013/3Loops*.tif')
spots_3p = []
# Loop through and find the spot number.
for f in images_3p:
    # Load the image.
    im = skimage.io.imread(f)

    # Feed it to our function.
    num_spots = spot_counter(im, threshold)

    #Store num_spots
    spots_3p.append(num_spots)


# Plot the spot number as function of time.
time_5p = np.arange(0, len(spots_5p), 1) * 10  # in seconds
time_3p = np.arange(0, len(spots_3p), 1) * 10  # in seconds

# Generate our plot.
plt.figure()
plt.plot(time_5p, spots_5p, label="5' labeled mRNA")
plt.plot(time_3p, spots_3p, label="3' labeled mRNA")
plt.xlabel('time (s)')
plt.ylabel('number of spots')
plt.show()


# Plot the synchronized spot numbers
# Looking at the nuclei, anaphase for 5' starts at frame 99
# anaphase starts for 3' at frame 125.
plt.figure()

# Index the time and spots such that they are synchronized.
plt.plot(time_5p[98:], spots_5p[98:], label="5' labeled")
plt.plot(time_3p[124:], spots_3p[124:], label="3' labeled")

# add labels and show.
plt.xlabel('time (s)')
plt.ylabel('number of spots')
plt.legend()
plt.show()
