#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  9 13:05:26 2021

3D image segmentation first test in python

original code from the image analysis lecture 
https://www.youtube.com/watch?v=lRtGqc5r6O0&list=PL5ESQNfM5lc7SAMstEu082ivW4BDMvd0U&index=33

and from 
https://scikit-image.org/docs/dev/auto_examples/applications/plot_3d_image_processing.html

@author: jingkui.wang
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from skimage import exposure, io, util
#from skimage.data import cells3d

from skimage.io import imread
import matplotlib.pyplot as plt


data = imread('../raw_imagines/0210217_nodrug_LDNSB_1_01.tif')

print("shape: {}".format(data.shape))
print("dtype: {}".format(data.dtype))
print("range: ({}, {})".format(data.min(), data.max()))

# Report spacing from microscope
original_spacing = np.array([0.2900000, 0.0650000, 0.0650000])


#print(image.shape)
#plt.imshow(image[1::1])
#plt.show()

# Let us try and visualize the (3D) image with io.imshow.
try:
    io.imshow(data, cmap="gray")
except TypeError as e:
    print(str(e))


xx = data

data = xx[:,:,:, 3]

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





