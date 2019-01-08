# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import time as tt
# For image processing
import skimage.io
import skimage.morphology
import skimage.segmentation
import scipy.ndimage

# In this script, we will determine the spring constant for a trapped
# plastic bead. This bead is taken from a stock with a mean diameter of one
# micron. The images we were used were captured at 50 ms intervals. Our
# general strategy will be to segment the bead from each individual image and
# identify the centroid. We'll then determine the mean squared displacement.

plt.close('all')
# Let's start with one image.
bead_ims = skimage.io.imread('data/optical_tweezer/trapped_bead_5.2x_4_MMStack_Pos0.ome.tif')
bead_im = bead_ims[0]
plt.figure()
plt.imshow(bead_im, cmap=plt.cm.Greys_r)
plt.show()

# We can see that the bead is light on a black background. Let's try
# some simple thresholding of this image. We'll convert this to a float first then look at the histogram.
im_float = (bead_im - bead_im.min()) /(bead_im.max() - bead_im.min())
im_blur = skimage.filters.gaussian(im_float, 3)

plt.figure()
plt.hist(im_blur.flatten(), bins=1000)
plt.xlabel('pixel value')
plt.ylabel('counts')
plt.show()


## Looking at the histogram, it seems like a threshold value of greater than 0.2 might be appropriate.
thresh = 0.6
im_thresh = (im_blur < thresh)
plt.figure()
plt.imshow(im_thresh, cmap = plt.cm.Greys_r)
plt.show()

# Check segmentation.
im_copy = np.copy(im_float)
bounds = skimage.segmentation.find_boundaries(im_thresh)
im_copy[bounds] = 1.0
merge = np.dstack((im_copy, im_float, im_float))
plt.figure()
plt.imshow(merge)
plt.show()
##
# That seems pretty good! Now to get the centroid, all we have to do is label our object and get the region props. Let's first remove small stuff that might be segmented in the background.
im_large = skimage.morphology.remove_small_objects(im_thresh)
im_border = skimage.segmentation.clear_border(im_large)
im_lab = skimage.measure.label(im_border)
# Now we'll label and extract the properties.
props = skimage.measure.regionprops(im_lab)
centroid = props[0].centroid

# Now let's also plot the centroid here.
plt.figure()
plt.imshow(im_float.T, cmap=plt.cm.viridis)
plt.plot(centroid[0], centroid[1], 'o', label='centroid')
plt.show()

# We had to transpose the image in the last plot in order for the axes to
# agree with eachother. Let's write this as a function and do it for each image

def find_centroid(im, threshold):
    """
    Extracts the centroid of a trapped bead.
    """
    # Make the image a float.
    im_float = (im - im.min()) / (im.max() - im.min())
    im_blur = skimage.filters.gaussian(im_float, 3)

    # Apply the threshold.
    im_thresh = im_blur < threshold



    # Get rid of small things.
    im_large = skimage.morphology.remove_small_objects(im_thresh)
    im_border = skimage.segmentation.clear_border(im_large)

    # Make sure only one object is segmented.
    im_lab, num_obj = skimage.measure.label(im_border, return_num=True)
    print(num_obj)
    if num_obj > 1:
        print('multiple objects found! Returning the centroid for largest object.')

    # Compute and return the centroid.
    props = skimage.measure.regionprops(im_lab)

    # Get the index of the largest area.
    areas = [prop.area for prop in props]
    ind = np.argmax(areas)
    centroid = props[ind].centroid
    return centroid


# Loop through each image. Save the x and y position.
centroid_x, centroid_y = [], []

# Get the file names.
length = len(bead_ims)
for i in range(length):
    print(i)
    # Load the image and process.
    im = bead_ims[i]
    x, y = find_centroid(im, thresh)
    centroid_x.append(x)
    centroid_y.append(y)


plt.figure()
plt.plot(centroid_x, centroid_y)
plt.xlabel('x position')
plt.ylabel('y position')
plt.show()

# We'll compute the mean position in both dimensions.
ip_dist = 0.042
x_micron = np.array(centroid_x) * ip_dist
y_micron = np.array(centroid_y) * ip_dist
mean_x = np.mean(x_micron)
mean_y = np.mean(y_micron)


# Now the mean squared displacement.
msd_x = np.mean((x_micron - mean_x)**2)
msd_y = np.mean((y_micron - mean_y)**2)

kBT = 4.1E-3  # in pN micrometers
k_x = kBT / (msd_x)
k_y = kBT / (msd_y)

print('Trapping force in x is ' + str(k_x) + ' pN•uM')
print('Trapping force in y is ' + str(k_y) + ' pN•uM')

fig, ax = plt.subplots(2, 1, sharex=True)
time = np.arange(0, length, 1)
ax[0].plot(time, centroid_x, '-')
ax[0].set_ylabel('x position')
ax[1].plot(time, centroid_y, '-')
ax[1].set_ylabel('y position')
ax[1].set_xlabel('time (frames)')
plt.show()


plt.figure()
for i, im in enumerate(bead_ims):
    im_float = (im - im.min()) / (im.max() - im.min())
    im_blur = skimage.filters.gaussian(im_float, 3)
    im_thresh = im_blur < thresh
    im_large = skimage.morphology.remove_small_objects(im_thresh)
    im_border = skimage.segmentation.clear_border(im_large)
    im_copy = np.copy(im_float)
    bounds = skimage.segmentation.find_boundaries(im_border)
    im_copy[bounds] = 1.0
    merge = np.dstack((im_copy, im_float, im_float))
    if i < 10:
        num = '00' + str(i)
    elif i < 100:
        num = '0' + str(i)
    else:
        num = str(i)

    with sns.axes_style('white'):
        plt.imshow(merge, interpolation='nearest')
        plt.plot(centroid_y[i], centroid_x[i], 's')
        plt.tight_layout()
        ax = plt.gca()
        ax.set_frame_on(False)
        plt.xticks([])
        plt.yticks([])
        plt.savefig('outpt/merge_' + num + '.tif')
        plt.clf()
