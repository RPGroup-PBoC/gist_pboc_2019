# Import the necessary modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')
# For image processing
import skimage.io

# In this script, we will compare theory with experimental data. In class, we
# derived "The French Flag" model for patterning in Drosophila embryos. We have
# provoded two images of developing embryos with different dosages of teh
# bicoid morphogen. To test our theory, we'll determine the positioning of the
# cephalic furrow in wild-type (dosage = 1)  and a mutant (dosage = 0.5) by
# using mouse clicks. Let's begin by loading in the wt embryo image.
wt_embryo = skimage.io.imread('data/french_flag/wild_type.tif')
plt.figure()
plt.imshow(wt_embryo)

# There are a few things to realize about this image. First, the embry os
# tildted and the exact position of the cephalic furrow is hard to pin point.
# We will have to use some trigonometry to  determine the position of the
# anterior-posterior axis and the relative positioning of the cephalic furrow.
# TO begin, we'll click on the position of the anterior, the posterior, and the
# cephalic furrow. Wherever you want to define these points is fine so long as
# you are consistent for each embry. To record the positions of the mouse
# clicks, we'll use the Python 'plt.ginput' command.
plt.figure()
plt.imshow(wt_embryo, cmap=plt.cm.Greys_r)
plt.show()
clicks = plt.ginput(2) # Will record two clicks.

# Now we'll extract the proper values from the click vector.
wt_ax = clicks[0][0]
wt_ay = clicks[0][1]
wt_px = clicks[1][0]
wt_py = clicks[1][1]

# Now that we have the position of the posterior and anterior poles, we'll
# find the position of the cephalic furrow.
plt.figure()
plt.imshow(wt_embryo, cmap=plt.cm.Greys_r)
plt.show()
furrow_click = plt.ginput(1)
wt_fx = furrow_click[0][0]
wt_fy = furrow_click[0][1]

# To deterime the position of the furrow, we'll find the length of the vector
# between teh anterior and posterior poles as well as the length between the
# cephalic furrow and the anterior pole.
wt_ap = np.sqrt((wt_py - wt_ay)**2 + (wt_px - wt_ax)**2)  # AP distance.
wt_fa = np.sqrt((wt_fy - wt_ay)**2 + (wt_fx - wt_ax)**2)  # Furrow position

# Using some trigonometry, we can find the length of the fragment of the
# anterior-posterior that intersects with the cephalic furrow. Alpha is an
# interior angle of teh triangle with hypotenouse AP. Beta is an interorior
# angel in the triangle formed with the hypotenouse CFA. Gamma is an interior
# angle of the triangel with CFA as the hypotenouse and a portion of AP as one
# side.
alpha = np.arctan((wt_py - wt_ay) / (wt_px - wt_ax))
beta = np.arctan((wt_fy - wt_ay) / (wt_fx - wt_ax))
gamma = beta - alpha

# With this information, we can compute the absolute length (in pixesl) of the
# fragement of AP betweenthe anterior and the cephalic furrow as well as the
# relative length.
wt_abs_cf = np.cos(gamma) * wt_fa
wt_rel_cf = wt_abs_cf / wt_ap

# As a sanity check, we know that the cephalic furrow in the wt embryo should be about 30% of the length of the embry. Let's see what we calculated.
print('The wt CF position is ' + str(wt_rel_cf))

# So now that we know how to calculate the position, let's do teh same procedure to the Bicoid mutant. We could simplify this by defining a function, but for now let's rewrite the same protocol.


# Read the mutant embry.
mut_embryo = skimage.io.imread('data/french_flag/bcd_mut.tif')

# Record the positions.
plt.figure()
plt.imshow(mut_embryo, cmap=plt.cm.Greys_r)
plt.show()
clicks = plt.ginput(2)
bcd_ax = clicks[0][0]
bcd_ay = clicks[0][1]
bcd_px = clicks[1][0]
bcd_py = clicks[1][1]

# Find the cephalic furrow position.
plt.figure()
plt.imshow(mut_embryo, cmap= plt.cm.Greys_r)
plt.show()
furrow_click = plt.ginput(1)

bcd_fx = furrow_click[0][0]
bcd_fy = furrow_click[0][1]

# Determine the position of the furrow.
bcd_ap = np.sqrt((bcd_py - bcd_ay)**2 + (bcd_px - bcd_ax)**2)  # AP distance.
bcd_fa = np.sqrt((bcd_fy - bcd_ay)**2 + (bcd_fx - bcd_ax)**2)  # Furrow position


# Compute the positioning for the mutant.
alpha = np.arctan((bcd_py - bcd_ay) / (bcd_px - bcd_ax))
beta = np.arctan((bcd_fy - bcd_ay) / (bcd_fx - bcd_ax))
gamma = beta - alpha

# Compute the length of the anterior posterior axis.
bcd_abs_cf = np.cos(gamma) * bcd_fa
bcd_rel_cf = bcd_abs_cf / bcd_ap


# Now we have the positioning of the cephalic furrow fo the WT and the bicoid
# mutant. Let's compare this data with the theory. We'll compute the curve we
# defined in class give our WT cephalic furrow position and plot our data
# points o the same set o axes. We'll start by defining some parameters.
lam = 0.2   # Decay length. Don't name this lambda!
xcf = wt_rel_cf
dosages = np.linspace(0.4, 2, 500)
new_xcf = xcf + np.log(dosages) * lam  # Our prediction.

# Plot it!
plt.figure()
plt.plot(dosages, new_xcf, '-', label='prediction')
plt.plot(1, wt_rel_cf, 'o', label='wt')   # the wild type.
plt.plot(0.5, bcd_rel_cf, 'o', label='mutant')
plt.xlabel('bicoid dosage (fractional)')
plt.ylabel('cf position from anterior (relative units)')
#plt.xscale('log')
plt.show()

# What do you say? Does the data agree with the theory?
