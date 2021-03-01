import numpy as np 
from scipy.io import loadmat
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from fast_poisson import fast_poisson
import cv2
import matplotlib.pyplot as plt
#from fast_poisson import poisson_reconstruct.
from calibration import calibration
import time
import pdb 

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
    diff_temp2[:,:,0] = (diff_temp2[:,:,0] - cali.zeropoint[0])/cali.lookscale[0]
    diff_temp2[:,:,1] = (diff_temp2[:,:,1] - cali.zeropoint[1])/cali.lookscale[1]
    diff_temp2[:,:,2] = (diff_temp2[:,:,2] - cali.zeropoint[2])/cali.lookscale[2]
    diff_temp3 = np.clip(diff_temp2,0,0.999)
    diff = (diff_temp3*cali.bin_num).astype(int)

    # pdb.set_trace()
    # _ = plt.hist(np.ndarray.flatten(r), bins='auto')
    # plt.show()
    # plt.figure()
    # plt.imshow((diff_temp1-np.min(diff_temp1)).astype(np.uint8))
    # plt.figure(0)
    # plt.imshow(ref_blur.astype(np.uint8))
    # plt.figure(1)
    # plt.imshow(test_img.astype(np.uint8))
    # plt.figure(2)
    # plt.imshow(((diff_temp2 - np.min(diff_temp2))/100).astype(np.uint8))
    # plt.show()
    grad_img = table[diff[:,:,0], diff[:,:,1],diff[:,:,2], :]
    return grad_img
    
    
def show_depth(depth, figure_num):
    fig = plt.figure(figure_num)
    ax = fig.gca(projection='3d')
    # ax.set_aspect('equal')
    X = np.arange(0, depth.shape[1], 1)*4.76/80
    Y = np.arange(0, depth.shape[0], 1)*4.76/80
    X, Y = np.meshgrid(X, Y)
    ax.plot_surface(X, Y, depth, cmap=cm.jet)
#    fig.colorbar(surf, shrink=0.5, aspect=5)
#    plt.show()

def contact_detection(raw_image, ref_blur,marker_mask, kernel):
    diff_img = np.max(np.abs(raw_image.astype(np.float32) - ref_blur),axis = 2)
    contact_mask = (diff_img> 25).astype(np.uint8)  #*(1-marker_mask)
    contact_mask = cv2.dilate(contact_mask, kernel, iterations=1)
    contact_mask = cv2.erode(contact_mask, kernel, iterations=1)
    return contact_mask
    
def marker_detection(raw_image_blur):
    m, n = raw_image_blur.shape[1], raw_image_blur.shape[0]
    raw_image_blur = cv2.pyrDown(raw_image_blur).astype(np.float32)
    ref_blur = cv2.GaussianBlur(raw_image_blur, (25, 25), 0)
    diff = ref_blur - raw_image_blur
    diff *= 16.0
    # cv2.imshow('blur2', blur.astype(np.uint8))
    # cv2.waitKey(1)
    diff[diff < 0.] = 0.
    diff[diff > 255.] = 255.

    # diff = cv2.GaussianBlur(diff, (5, 5), 0)
    # cv2.imshow('diff', diff.astype(np.uint8))
    # cv2.waitKey(1)
    # mask = (diff[:, :, 0] > 25) & (diff[:, :, 2] > 25) & (diff[:, :, 1] >
    #                                                       120)
    mask_b = diff[:, :, 0] > 150 
    mask_g = diff[:, :, 1] > 150 
    mask_r = diff[:, :, 2] > 150 
    mask = (mask_b*mask_g + mask_b*mask_r + mask_g*mask_r)>0
    # cv2.imshow('mask', mask.astype(np.uint8) * 255)
    # cv2.waitKey(1)
    mask = cv2.resize(mask.astype(np.uint8), (m, n))
    return mask 

def make_kernal(n,k_type):
    if k_type == 'circle':
        kernal = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(n,n))
    else:
        kernal = cv2.getStructuringElement(cv2.MORPH_RECT,(n,n))
    return kernal 
#
    
if __name__ == '__main__': 
    table2 = np.load('table_smooth.npy')
    abe_array = np.load('abe_corr.npz')
    x_index, y_index = abe_array['x'], abe_array['y']
    kernel1 = make_kernal(3,'circle')
    kernel2 = make_kernal(25,'circle')

    cali = calibration()
    pad = 20
    ref_img = cv2.imread('./test_data/ref.jpg')
    ref_img = ref_img[x_index, y_index, :]
    test_img = cv2.imread('./test_data/gelslim1_250.jpg')
    test_img = test_img[x_index, y_index, :]
#    ref_img = test_img.copy()
    ref_img = cali.crop_image(ref_img, pad) 
    marker = cali.mask_marker(ref_img)
    keypoints = cali.find_dots(marker) 
    marker_mask = cali.make_mask(ref_img, keypoints)
    marker_image = np.dstack((marker_mask, marker_mask, marker_mask))
    ref_img = cv2.inpaint(ref_img,marker_mask,3,cv2.INPAINT_TELEA)
    red_mask = (ref_img[:,:,2] > 12).astype(np.uint8)
    ref_blur = cv2.GaussianBlur(ref_img.astype(np.float32), (3, 3), 0) + 1
    # pdb.set_trace()
#    ref_blur_small = cv2.pyrDown(ref_blur).astype(np.float32)
    blur_inverse = 1 + ((np.mean(ref_blur)/ref_blur)-1)*2;
    test_img = cali.crop_image(test_img, pad)
    test_img = cv2.GaussianBlur(test_img.astype(np.float32), (3, 3), 0)
#    t1 = time.time()
    marker_mask = marker_detection(test_img)
    marker_mask = cv2.dilate(marker_mask, kernel1, iterations=1)
    contact_mask = contact_detection(test_img, ref_blur,marker_mask, kernel2)

    # mask_2_show = np.dstack((np.zeros_like(marker_mask), marker_mask, np.zeros_like(marker_mask)))*40 +  test_img.astype(np.uint8)
    # plt.figure(20)
    # plt.imshow(mask_2_show)
    # plt.figure(21)
    # plt.imshow(contact_mask)
#    plt.show()
    
    grad_img2 = matching_v2(test_img, ref_blur, cali, table2, blur_inverse)

    grad_img2[:,:,0] = grad_img2[:,:,0] * (1-marker_mask) * red_mask
    grad_img2[:,:,1] = grad_img2[:,:,1] * (1-marker_mask) * red_mask

#    depthgx = np.array(gx.ImGradX)1 = fast_poisson(grad_img1[:,:,0], grad_img1[:,:,1])
    depth2 = fast_poisson(grad_img2[:,:,0], grad_img2[:,:,1])
   # depth1[depth1<0] = 0
    depth2[depth2<0] = 0
   # show_depth(depth1,99)
    show_depth(depth2,100)
   # print(time.time()-t1)
    plt.figure(0)
    plt.imshow(grad_img2[:,:,0])
    plt.figure(1)
    plt.imshow(grad_img2[:,:,1])
    plt.figure(2)
    plt.imshow(depth2)
   # plt.figure(3)
   # plt.imshow((ref_blur)/255.)
    # plt.figure(5)
    # plt.imshow(cv2.cvtColor((test_img-ref_blur)/70, cv2.COLOR_BGR2RGB))
    plt.show()
   # cv2.imshow('diff',(((test_img-ref_blur)+150)/400*255).astype(np.uint8))
   # cv2.waitKey(0)
#%%
# cv2.imshow('test_image', test_img.astype(np.uint8))
# cv2.waitKey()

    