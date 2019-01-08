# Import the important modules.
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import utilities for image processing.
import skimage.io
sns.set_style('white')

# This script will introduce you to the syntax and various operations in the
# Python programming language as well as introduce you to some basic image
# processing. We'll learn this by 'sizing up' some E. coli cells. You should
# note that there are many different programming languages (C, C++, Matlab,
# python, java, julia, etc) that are very useful for scientists of every
# discipline. While we will use Python for this course, I urge you to not
# subscribe to any programming language with religious fervor.

# There is no standard GUI for Python. In this course, we will use a
# combination of the Atom text editor and the IPython interpreter console.
# Using the Anaconda python distribution, we can launch the latter from the
# 'Navigator' application. In the IPython interpreter (also called the
# qtconsole), we can run small snippets of code, individual lines, and even
# entire Python scripts. To get a flavor for the Syntax and behavior of Python,
# Type the following commands into the IPython interpreter and see what they
# do. We will make these 'print' in this script so the output will be generated
# each time we run it.

# Mathematical operations behave almost exactly as you would expect.
print(1 + 1)  # This should give us 2
print(2 * 8)  # This should give us 16
print(2**2)  # This should give us 4.

# Note that exponentiation in Python is ** not ^.
# We can access other mathematical operations using the Numerical Python
# (NumPy) module which we have imported at the beginning of this script. Note
# that we imported it as `np`. This means that we can access all functions and
# methods in the numpy package by simply prefacing it with np. For example,

print(np.sin(2 * np.pi))  # This should give us 2.
print(np.exp(2 * np.abs(-5)))  # This should give us e^10
print(np.log(np.exp(4)))  # This should give us 4.
print(np.log10(1E4))  # This will take the log base 10 and give 4.

# We can also store values as variables. We assign these values by using an
# equals sign (=).
a = 2
b = 3
c = a * b / a
print(c)  # Should give 3.

# We can store sets of values as lists or arrays.
list_values = [1, 2, 3, 4]
array_values = np.array([1, 2, 3, 4])

# We can index and extract individual values from this lsit and array.
print(list_values[0])  # This should regurn 1.
print(array_values[2])  # This should return 3.
print(array_values[-1])  # This will return the last value of the array.

# Note that in Python, indexing begins at zero.
# Using numpy, we can also automatically generate specific kinds of arrays.
linear_array = np.linspace(0, 10, 100)
print(linear_array)  # 100 points starting at 0 and ending 10.
log_array = np.logspace(-2, -1, 100)
print(log_array)  # 100 points starting at 1E-2 ending 1E-4
arranged_array = np.arange(0, 10, 1)
print(arranged_array)  # Values between 0 and 10 taking steps of 1.

# Now that we have some of the basic syntax down, let's begin by reading in an
# image of a graticule, which is basically a ruler used for microscopy. To
# read in the image, we will use another module called skimage.io.

# Remember, an image is just data -- a simple two-dimensional array in which
# each element corresponds to a pixel value.

# Read in the image.
grat_im = skimage.io.imread('data/sizing_up_ecoli/Graticule100x.tif')

# If we just print out the values of the grat_im variable, we will see that is
# just a simple array of values.
print(grat_im)

# This isn't really a useful way to display an image, however. Let's show
# the image in the plotting window. To do this, we will use yet another module
# called matplotlib.pyplot that we imported as plt.

# Show the graticule image.
plt.figure()
plt.imshow(grat_im)
plt.show()

# Each major division on this image is 10 microns apart. If we consider this
# image as an array, we see that if we were to go across a given row of the
# image, therew ould be a periodic osciallation inthe pixel values. Dark values
# would correspond to the ticks in the graticules while the white spaces would
# correspond to the spaces. Using this, let's extract the pixel value
# information for a single row in the image. Let's figure out what the size
# of the image is.
print(np.shape(grat_im))


# This tells us that the image is 982 rows by 1311 coluns. Let's look at the
# pixel values for row 491, exactly down the middle.
line_profile = grat_im[491, :]

# We'll plot these values agains thte length of the xaxis of the image.
x = np.arange(0, len(line_profile), 1)

# Now well plot it as another figure.
plt.figure()
plt.plot(x, line_profile, 'r-')  # plotting the profile as a red line.
plt.xlabel('x position')
plt.ylabel('pixel value')

# Just by zooming in on the image, we can see that there is a valley in the
# pixel value at about pixel number 100 and another at 250. This means that
# by our crude estimate, the distance between pixels is as follows.
crude_ip_dist = 10 / 150
print(crude_ip_dist)

# which is about 66nm per pixel. We can be a little more creative and actually
# click on the positions which we think are the middle of the two major
# divisions.

# Close all the other figures and plot the graticule image.
plt.close('all')
plt.figure()
plt.imshow(grat_im)
clicks = plt.ginput(2)  # The 2 tells it will wait for two clicks.

# These clicks store the x and y position for each pixel that we clicked on.
# We can easily compute the interpixel distance by taking the difference
# getween the two x positions, since y position doesn't matter.
diff_clicks = clicks[1][0] - clicks[0][0]
clicks_ip_dist = 10 / diff_clicks
print(click_ip_dist)

# This interpixel distance is about the same as we got by our crude measure
# made above. Now that we have an idea of what the distance between pixels is,
# let's make a measurement of the length of these cells.

# Load int the image of the cells.
cell_im = skimage.io.imread('data/sizing_up_ecoli/Ecoli100x.tif')

# Let's take a look at it for good measure.
plt.figure()
plt.imshow(cell_im)
cell_clicks = plt.ginput(2)

# We can now determine the length of the bacterium using the pythagorean
# theorem.
cell_length_pix = np.sqrt((cell_clicks[0][0] - cell_clicks[1][0])**2 + \
                          (cell_clicks[0][1] - cell_clicks[1][1])**2)

# Now it is a simple multiplication to get the physical distance.
cell_length_micron = clicks_ip_dist * cell_length_pix
print('Our measured cell is ' + str(cell_length_micron))

# That's not too bad! Let's upt a 10 micron scale bar on our image.
# We can do that by changing a row of pixels 10 microns long on the image
# all to balck.
scale_bar_length = 10 / clicks_ip_dist
copy_im = np.copy(cell_im)  # Make a copy for safety.

# Change the row of pixels to black.
copy_im[900:920, 20:20 + np.round(scale_bar_length)] = 0

# Show the image.
plt.figure()
plt.imshow(copy_im)

# Let's save it for fun.
plt.savefig('data/sizing_up_ecoli/ecoli_with_scalebar.tif')
