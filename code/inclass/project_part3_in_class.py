# Import our awesome modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn
import glob

# Image processing modules.
import skimage.io
import skimage.filters
import skimage.measure
import skimage.segmentation

# Here's what we've done so far.
im = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')
yfp_im = skimage.io.imread('data/lacI_titration/O2_delta_yfp_pos_16.tif')

# Normalize the image.
im_norm = (im - im.min()) / (im.max() - im.min())

# Do the background subtraction
im_blur = skimage.filters.gaussian(im_norm, 50.0)
im_sub = im_norm - im_blur

# Threshold the image.
im_thresh = im_sub < -0.2

# Label our image.
im_label = skimage.measure.label(im_thresh)
props = skimage.measure.regionprops(im_label)

# We want to keep the cells with a given area.
approved_objects = np.zeros_like(im_label)
ip_dist = 0.160  # in units of microns per pixel
for prop in props:
    obj_area = prop.area * ip_dist**2
    if (obj_area > 0.5) & (obj_area < 5):
        approved_objects += (im_label == prop.label)

# Extract the intensities.
mean_int = []
im_relab = skimage.measure.label(approved_objects)
props = skimage.measure.regionprops(im_relab, intensity_image=yfp_im)
for prop in props:
    mean_int.append(prop.mean_intensity)


plt.figure()
plt.hist(mean_int, bins=10)
plt.xlabel('mean pixel intensity')
plt.ylabel('count')
plt.show()


def phase_segmentation(image, threshold):
    """
    Performs segmentation on a phase image.
    """
    # Normalize the image
    im_norm = (image - image.min()) / (image.max() - image.min())

    # Do a background subtraction
    im_blur = skimage.filters.gaussian(image, 50.0)
    im_sub = im_norm - im_blur

    # Threshold the image
    im_thresh = im_sub < -0.2

    # Label the image
    im_label = skimage.measure.label(im_thresh)

    # Get the properties and apply an area threshold
    props = skimage.measure.regionprops(im_label)

    # Make an empty image to store the approved cells
    approved_objects = np.zeros_like(im_label)

    # Apply the area filters
    for prop in props:
        obj_area = prop.area * 0.160**2  # Given the interpixel distance
        if (obj_area > 0.5) & (obj_area < 5):
            approved_objects += (im_label==prop.label)

    # Relabel the image.
    return im_relab

def extract_intensity(mask, yfp_image):
    """
    Extract the mean intensity from a segmented image.
    """
    # Get the region properties for the image.
    props = skimage.measure.regionprops(mask, intensity_image=yfp_image)

    # Make a vector to store the mean intensities
    mean_int = []
    for prop in props:
        intensity = prop.mean_intensity
        mean_int.append(intensity)

    return mean_int


# With these functions in hand, let's loop over autofluorescence and delta.
delta_phase = glob.glob('data/lacI_titration/O2_delta_phase*.tif')
delta_yfp = glob.glob('data/lacI_titration/O2_delta_yfp_pos*.tif')

delta_mean_int = []
for i in range(len(delta_phase)):
    im = skimage.io.imread(delta_phase[i])
    yfp_im = skimage.io.imread(delta_yfp[i])

    # Put it through our functions.
    mask = phase_segmentation(im, -0.2)
    ints = extract_intensity(mask, yfp_im)

    #Loop through the intensity and add it.
    for value in ints:
        delta_mean_int.append(value)

# Now do the same for the autoflurescent samples.
auto_phase = glob.glob('data/lacI_titration/O2_auto_phase_*.tif')
auto_yfp = glob.glob('data/lacI_titration/O2_auto_yfp_*.tif')
auto_mean_int = []
for i in range(len(auto_phase)):
    im = skimage.io.imread(auto_phase[i])
    yfp_im = skimage.io.imread(auto_yfp[i])
    mask = phase_segmentation(im, -0.2)
    ints = extract_intensity(mask, yfp_im)
    for value in ints:
        auto_mean_int.append(value)

# Now generate the histograms of each.
plt.figure()
plt.hist(delta_mean_int, bins=100)
plt.xlabel('mean pixel intensity')
plt.ylabel('counts')
plt.title('delta sample')

plt.figure()
plt.hist(auto_mean_int, bins=100)
plt.xlabel('mean pixel intensity')
plt.ylabel('counts')
plt.title('autofluorescent sample')
plt.show()
