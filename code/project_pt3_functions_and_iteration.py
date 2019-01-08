# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Image processing utilities.
import skimage.io
import skimage.measure
import skimage.filters
import skimage.segmentation

# For file management.
import glob

# In this script, we'll write all of the image processing steps we've executed
# so far as two functions. We would like to split these functions up in to two
# parts.
#  1. Performing the segmentation. This should take an image, a threshold
#    value, and some bounds for the area filtering. This will return the
#    labeled segmentation  mask.
#  2. Extracting the fluorescence intentsities of each object. This should take
#    a segmentation mask and a fluorescence image. This will return an array of
#    the mean fluorescence for each cell.

# Why would we want to write this as a function? We'll be doing this same
# procedure many times over each image in our sample. We'll have more than 100
# images, so it would be best to have this be modular to save our fingers some
# work. To verify that our functions work, we'll make sure that we get the
# exact same result when type each command out by hand on a single image.


# To start, let's retype what we've done so far. We'll generate a final
# segmentation mask and a list of cell intensities.
phase_image = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')

# Normalize the image and perform a background subtraction.
im_float = (phase_image - phase_image.min()) / (phase_image.max() - phase_image.min())
im_blur = skimage.filters.gaussian(im_float, sigma=50.0)
im_sub = im_float - im_blur

# Apply the threshold and the area filter.
im_thresh = im_sub < -0.2
im_label = skimage.measure.label(im_thresh)
props = skimage.measure.regionprops(im_label)
ip_dist = 0.160  # in units of microns per pixel.
approved_objects = np.zeros_like(im_label)
for prop in props:
    area = prop.area * ip_dist**2
    if (area > 0.5) & (area < 6):
        approved_objects += im_label==prop.label

# Clear the border and relabel.
im_border = skimage.segmentation.clear_border(approved_objects > 0)
im_relab = skimage.measure.label(im_border)

# Show the final segmentation mask.
plt.figure()
plt.imshow(im_relab, cmap=plt.cm.spectral_r)
plt.show()

# Load the fluorescence image and compute the cell intensities.
fluo_im = skimage.io.imread('data/lacI_titration/O2_delta_yfp_pos_16.tif')
props = skimage.measure.regionprops(im_relab, intensity_image=fluo_im)
cell_intensities = []
for prop in props:
    cell_intensities.append(prop.mean_intensity)

# Print the mean cell intensities to the screen.
mean_int = np.mean(cell_intensities)
print("The mean cell intensity coded by hand is " + str(mean_int) + " counts.")

# This looks like what we've done in the past three project parts. Let's try
# writing these as functions. We'll write this functions such that they take
# only one or two arguments with the other parameters as keyword arguments.


# We'll start with the segmentation function.
def phase_segmentation(image, threshold, area_bounds=[0.5, 6.0],
                       ip_dist=0.160):
    """
    Segement a phase image and return the mask.

    Parameters
    ----------
    image : 2d-array
        The phase image to be segmented. This image will be converted to a
        float type.
    threshold : float
        Threshold value for the segmentation. This function will select objects
        below this threshold value.
    area_bounds : list, default=[0.5, 6.0]
        Area bounds for identified objects. This should be a list of two entries.
    ip_dist : int or float, default = 0.160
        Interpixel distance for the camera. This should be in units of microns
        per pixel.

    Returns
    -------
    final_seg : 2d-array
        Final, labeled segmentation mask.
    """

    # First is to convert the image to a float.
    im_float = (image - image.min()) / (image.max() - image.min())

    # Do a background subtraction.
    im_blur = skimage.filters.gaussian(im_float, sigma=50.0)
    im_sub = im_float - im_blur

    # Apply the threshold.
    im_thresh = im_sub < threshold  # Note that we are using the provided arg

    # Label the image and apply the area bounds.
    im_lab = skimage.measure.label(im_thresh)
    props = skimage.measure.regionprops(im_lab)
    approved_objects = np.zeros_like(im_lab)
    for prop in props:
        area = prop.area * ip_dist**2
        if (area > area_bounds[0]) & (area < area_bounds[1]):
            approved_objects += im_lab == prop.label

    # Clear the border and relabel.
    im_border = skimage.segmentation.clear_border(approved_objects > 0)
    final_seg = skimage.measure.label(im_border)

    # Return the final segmentation mask
    return final_seg


# Now let's try writing one for to extract the mean intensities.
def extract_intensity(seg, fluo_im):
    """
    Extracts the mean intensity of objects in a segmented image.

    Parameters
    ----------
    seg : 2d-array, int
        Segmentation mask with labeled objects.
    fluo_im : 2d-array, int
        Fluorescence image to extract intensities from.

    Returns
    -------
    cell_ints : 1d-array
        Vector of mean cell intensities. This has a length the same as the
        number of objects in the provided segmentation mask.
    """

    # Get the region props of the fluorescence image using the segmentation
    # mask.
    props = skimage.measure.regionprops(seg, intensity_image=fluo_im)
    cell_ints = []
    for prop in props:
        cell_ints.append(prop.mean_intensity)

    # Convert the cell_ints to an array and return.
    return np.array(cell_ints)


# Let's test these two functions out on our image and make sure that they
# return the same values that we got from our hand coded version.
seg_mask = phase_segmentation(phase_image, -0.2)
cell_ints = extract_intensity(seg_mask, fluo_im)
mean_int_func = np.mean(cell_ints)
print("The mean cell intensity from our function is " + str(mean_int_func) + " counts.")


# Finally, we can test that the segmentation masks are the same by checking
# if each pixel in the two masks are identical.
print((seg_mask == im_relab).all())  # Checks that all pixels are the same.


# Great! They all seem to work. Ultimately, we would like to iterate this over
# all of our images for each concentration and each channel. Let's try doing
# this with our functions and generate a distribution of the fluorescence
# intensity for all images in the O2 delta images. To get a list of all of
# the images in the directory, we will use the glob module.
phase_names = glob.glob('data/lacI_titration/O2_delta_phase*.tif')
fluo_names = glob.glob('data/lacI_titration/O2_delta_yfp*.tif')

# We'll make an empty vector to store all of our cell intensity vectors.
cell_intensities = []

# Now we just have to loop through each file.
for i in range(len(phase_names)):
    # Load the two images.
    phase_im = skimage.io.imread(phase_names[i])
    fluo_im = skimage.io.imread(fluo_names[i])

    # Perform the segmentation.
    seg_mask = phase_segmentation(phase_im, -0.2)

    # Extract the intensities and append them.
    ints = extract_intensity(seg_mask, fluo_im)
    for val in ints:
        cell_intensities.append(val)


# Now we can plot the distribution!
print('Extracted ' + str(len(cell_intensities)) + ' intensity values.')
plt.figure()
plt.hist(cell_intensities, bins=100)
plt.xlabel('mean YFP cell intensity (a. u.)')
plt.ylabel('counts')
plt.title('constitutive expression')
plt.show()

# We see that we quite a lot of cells. Their distribution has a range between
# 2000 and 6000, like we saw in the previous distribution from a single image.

# Let's look at another image set to see how well our functions work. Let's
# take a look at the fluorescence distribution for the O2 autofluorescent
# strain. Remember, these cells should be very dark compared to the delta
# strain.
phase_names_auto = glob.glob('data/lacI_titration/O2_auto_phase*.tif')
fluo_names_auto = glob.glob('data/lacI_titration/O2_auto_yfp*.tif')
auto_cell_intensities = []
for i in range(len(phase_names_auto)):
    phase_im = skimage.io.imread(phase_names_auto[i])
    fluo_im = skimage.io.imread(fluo_names_auto[i])
    seg_mask = phase_segmentation(phase_im, -0.2)
    ints = extract_intensity(seg_mask, fluo_im)
    for val in ints:
        auto_cell_intensities.append(val)

# Let's look at the distribution.
# Now we can plot the distribution!
print('Extracted ' + str(len(auto_cell_intensities)) + ' auto intensity values.')
plt.figure()
plt.hist(auto_cell_intensities, bins=100)
plt.xlabel('mean YFP cell intensity (a. u.)')
plt.ylabel('counts')
plt.title('autofluorescence')
plt.show()

# That's quite a difference!

# In our next script, we will iterate this over all of the images in our
# collection and test our theory for fold-change in gene expression for simple
# repression.
