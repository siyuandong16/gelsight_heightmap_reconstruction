import numpy as np 
#from scipy.io import loadmat
from mpl_toolkits.mplot3d import Axes3D
from fast_poisson import fast_poisson
import cv2
import matplotlib.pyplot as plt
#from fast_poisson import poisson_reconstruct.
from calibration import image_processor, calibration
import time

def matching(test_img, ref_blur,cali,table):
    diff = test_img - ref_blur
    
    diff[:,:,0] = np.clip((diff[:,:,0] - cali.blue_range[0])*cali.ratio, 0, cali.blue_bin-1)
    diff[:,:,1] = np.clip((diff[:,:,1] - cali.green_range[0])*cali.ratio, 0, cali.green_bin-1)
    diff[:,:,2] = np.clip((diff[:,:,2] - cali.red_range[0])*cali.ratio, 0, cali.green_bin-1)
    diff = diff.astype(int)
    grad_img = table[diff[:,:,0], diff[:,:,1],diff[:,:,2], :]
    return grad_img

def matching_v2(test_img, ref_blur,cali,table, blur_inverse):
    
    diff_temp1 = test_img - ref_blur
    diff_temp2 = diff_temp1 * blur_inverse
    diff_temp3 = np.clip((diff_temp2-cali.zeropoint)/cali.lookscale,0,0.999)
    diff = (diff_temp3*cali.bin_num).astype(int)
    grad_img = table[diff[:,:,0], diff[:,:,1],diff[:,:,2], :]
    return grad_img
    
    
    

    
if __name__ == '__main__': 
    table = np.load('table_2.npy')
#    plt.imshow(table[:,:,80,0])
#    plt.show()
    imp = image_processor()
    cali = calibration()
    pad = 20
    ref_img = cv2.imread('./test_data/ref.jpg')
    test_img = cv2.imread('./test_data/sample2_1010.jpg')
#    ref_img = test_img.copy()
    ref_img = imp.crop_image(ref_img, pad)
    ref_blur = cv2.GaussianBlur(ref_img.astype(np.float32), (49, 49), 49)
    blur_inverse = 1+((np.mean(ref_blur)/ref_blur)-1)*2;
    test_img = imp.crop_image(test_img, pad)
    test_img = cv2.GaussianBlur(test_img.astype(np.float32), (3, 3), 0)
    t1 = time.time()
    grad_img = matching_v2(test_img, ref_blur, cali, table, blur_inverse)
#    grad_img = matching(test_img, ref_blur, cali, table)
    depth = fast_poisson(grad_img[:,:,0], grad_img[:,:,1])
    print(time.time()-t1)
    plt.figure(0)
    plt.imshow(grad_img[:,:,0])
    plt.figure(1)
    plt.imshow(grad_img[:,:,1])
    plt.figure(2)
    plt.imshow(depth)
    plt.figure(3)
    plt.imshow((ref_blur)/255.)
    plt.figure(4)
    plt.imshow((test_img-ref_blur)/70.)
    plt.show()
    


#%%

    