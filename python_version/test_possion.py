import numpy as np 
from scipy.io import loadmat
from fast_possion import fast_poisson
import cv2
import matplotlib.pyplot as plt
from fast_poisson import poisson_reconstruct

gradx = loadmat('im_gradx.mat')
grady = loadmat('im_grady.mat')


imgradx = np.array(gradx['ImGradX'])
imgrady = np.array(grady['ImGradY'])

m,n = imgradx.shape
height_map = fast_poisson(imgradx, imgrady)

height_map2 = poisson_reconstruct(imgrady, imgradx, np.zeros(imgradx.shape))


plt.imshow(height_map)
plt.figure(2)
plt.imshow(height_map2)
plt.show()
# print(gradx.keys())
# print(grady.keys())
# print(imgradx.shape)
# print(imgradx[50,50])