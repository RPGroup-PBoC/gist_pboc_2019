# Import our stuff.
import numpy as np
import matplotlib.pyplot as plt
import seaborn

# For image processing
import skimage.io
import skimage.filters
import skimage.segmentation
import skimage.measure

# Do the same procedure as yesterday.
image = skimage.io.imread('data/lacI_titration/O2_delta_phase_pos_16.tif')

# Renormalzie the image.
im_norm = (image - image.min()) / (image.max() - image.min())

# Perform a background subtraction.
im_blur = skimage.filters.gaussian(im_norm, 50.0)
im_sub = im_norm - im_blur

# Threshold our image using yesterday's value.
im_thresh = im_sub < -0.2
plt.figure()
plt.imshow(im_thresh, cmap='gray')
plt.show()

# Label the image and print the number of segmented objects.
im_lab, num_objects = skimage.measure.label(im_thresh, return_num=True)
print('We found ' + str(num_objects) + ' objects. Is this right?')

# Show the labeled image.
plt.figure()
plt.imshow(im_lab, cmap='spectral')
plt.show()

# Extract the areas and plot a histogram.
plt.close('all')
props = skimage.measure.regionprops(im_lab)

# Get the areas from regionprops.
ip_dist = 0.160  # physical distance of pixels in microns per pixel
areas = []
for prop in props:
    obj_area = prop.area * ip_dist**2
    areas.append(obj_area)

plt.figure()
plt.hist(areas, bins=75)
plt.xlabel('area in pixels')
plt.ylabel('counts')
plt.show()

# Keep objects between 1 micron squared and 5 microns squared
approved_objects = np.zeros_like(im_lab)

# Iterate through each property and assess the area
for prop in props:
    obj_area = prop.area * ip_dist**2
    if (obj_area > 1) & (obj_area < 5):
        approved_objects = approved_objects + (im_lab==prop.label)

# Relabel our cells
approved_objects = approved_objects > 0
im_relab, new_count= skimage.measure.label(approved_objects, return_num=True)
print('Our new count is ' + str(new_count) + ' cells. Hooray!')

# Show the relabeled image.
plt.figure()
plt.imshow(im_relab, cmap='spectral')
plt.show()

# Get the fluorescence information of these cells.
yfp_im = skimage.io.imread('data/lacI_titration/O2_delta_yfp_pos_16.tif')
yfp_props = skimage.measure.regionprops(im_relab, intensity_image=yfp_im)
mean_intensities = []
for prop in yfp_props:
    mean_int = prop.mean_intensity
    mean_intensities.append(mean_int)

# Show the histogram of fluorescence
plt.figure()
plt.hist(mean_intensities, bins=10)
plt.xlabel('mean pixel intensity (a. u.)')
plt.ylabel('counts')
plt.show()
