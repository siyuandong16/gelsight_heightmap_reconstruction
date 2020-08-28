import cv2
import numpy as np 
import matplotlib.pyplot as plt  
import glob
from fast_poisson import fast_poisson
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
from skimage.restoration import inpaint




table = np.load('table_3.npy')
table_smooth = np.load('table_3_smooth.npy') 
count_map = np.load('count_map_3.npy')

print('table size', table.shape)
print('countmap size', count_map.shape)

plt.figure(0)
plt.imshow(table[:,:,45,0])

plt.figure(1)
plt.imshow(table_smooth[:,:,45,0])

plt.figure(2)
plt.imshow(count_map[:,:,45])
plt.show()
