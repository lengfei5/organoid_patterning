#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  9 13:05:26 2021

3D image segmentation first test in python

original code from the image analysis lecture 
https://www.youtube.com/watch?v=lRtGqc5r6O0&list=PL5ESQNfM5lc7SAMstEu082ivW4BDMvd0U&index=33

and from 
https://scikit-image.org/docs/dev/auto_examples/applications/plot_3d_image_processing.html

and also from 
https://github.com/scikit-image/skimage-tutorials/blob/main/lectures/solutions/4_segmentation.ipynb

@author: jingkui.wang
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from skimage import exposure, io, util
#from skimage.data import cells3d

from skimage.io import imread
import matplotlib.pyplot as plt

import skimage.segmentation as seg
from skimage import filters
from skimage import draw
from skimage import color
from skimage import exposure


import cv2
import os
import glob

#import itk
#from itkwidgets import view

# test itk 
#image = itk.imread('../raw_imagines/0210217_nodrug_LDNSB_1_01.tif')
#print(image)
#view(image)

image = imread('../raw_imagines/0210217_nodrug_LDNSB_1_01.tif')
data = image

print("shape: {}".format(data.shape))
print("dtype: {}".format(data.dtype))
print("range: ({}, {})".format(data.min(), data.max()))

# Report spacing from microscope
original_spacing = np.array_type([0.2900000, 0.0650000, 0.0650000])


#print(image.shape)
#plt.imshow(image[1::1])
#plt.show()

# Let us try and visualize the (3D) image with io.imshow.
try:
    io.imshow(data, cmap="gray")
except TypeError as e:
    print(str(e))




data = image[:,:,:, 0]

def show_plane(ax, plane, cmap="viridis", title=None):
    ax.imshow(plane, cmap=cmap)
    ax.axis("off")

    if title:
        ax.set_title(title)


(n_plane, n_row, n_col) = data.shape

_, (a, b, c) = plt.subplots(nrows = 1, ncols=3, figsize=(15, 5))

show_plane(a, data[n_plane // 2], title=f'Plane = {n_plane // 2}')
show_plane(b, data[10,:,:], title=f'Row = {n_row // 2}')
show_plane(c, data[15], title=f'Column = {n_col // 2}')


# As hinted before, a three-dimensional image can be viewed as a series of two-dimensional planes. 
# Let us write a helper function, display, to display 30 planes of our data. By default, every other plane is displayed.
def display(im3d, cmap="gray", step=5):
    _, axes = plt.subplots(nrows=4, ncols=5, figsize=(16, 14))

    vmin = im3d.min()
    vmax = im3d.max()

    for ax, image in zip(axes.flatten(), im3d[::step]):
        ax.imshow(image, cmap=cmap, vmin=vmin, vmax=vmax/5)
        ax.set_xticks([])
        ax.set_yticks([])

display(image[:,:,:, 0])
display(image[:,:,:, 1])
display(image[:,:,:, 2])
display(image[:,:,:, 3])

def plot_hist(ax, data, title=None):
    # Helper function for plotting histograms
    ax.hist(data.ravel(), bins=256)
    ax.ticklabel_format(axis="y", style="scientific", scilimits=(0, 0))

    if title:
        ax.set_title(title)

data = image[:,:,:, 0]
gamma_low_val = 0.1
gamma_low = exposure.adjust_gamma(data, gamma=gamma_low_val)

gamma_high_val = 0.3
gamma_high = exposure.adjust_gamma(data, gamma=gamma_high_val)

_, ((a, b, c), (d, e, f)) = plt.subplots(nrows=2, ncols=3, figsize=(12, 8))


show_plane(a, data[32], title='Original')
show_plane(b, gamma_low[32], title=f'Gamma = {gamma_low_val}')
show_plane(c, gamma_high[32], title=f'Gamma = {gamma_high_val}')

plot_hist(d, data)
plot_hist(e, gamma_low)
plot_hist(f, gamma_high)

# test some segmentation methods
def image_show(image, nrows=1, ncols=1, cmap='gray', **kwargs):
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(16, 16))
    ax.imshow(image, cmap='gray')
    ax.axis('off')
    return fig, ax

data = image[20:,:,:,3]
fig, ax = plt.subplots(1, 1)
ax.hist(data.ravel(), bins=500)
ax.set_xlim(0, 1000);

# supervised threshod
data = image[:,:,:,3]
data_segmented = data > 300 # Your code here

display(data_segmented)

#image_show(data_segmented);

# test watershed algorithm in 3D 
from skimage.segmentation import watershed
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
# Now we want to separate the two objects in image
# Generate the markers as local maxima of the distance to the background
distance = ndi.distance_transform_edt(data_segmented)

coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=data_segmented)
mask = np.zeros(distance.shape, dtype=bool)
mask[tuple(coords.T)] = True
markers, _ = ndi.label(mask)
labels = watershed(-distance, markers, mask=image)


# Unsupervised segmentation
#data_slic = seg.slic(data, n_segments=300)

#image_show(color.label2rgb(data_slic, data, kind='overlay'));

#display(data_slic)

#image_show(data_slic)

# test chan_verse method (only for 2D)
#chan_vese = seg.chan_vese(data)

#fig, ax = image_show(data)
#ax.imshow(chan_vese == 0, alpha=0.3);

