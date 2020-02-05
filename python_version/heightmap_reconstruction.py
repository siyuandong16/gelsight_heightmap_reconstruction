#!/usr/bin/env python


from sensor_msgs.msg import CompressedImage
from std_msgs.msg import Bool
import numpy as np
import time
from scipy import ndimage
import matplotlib.pyplot as plt
from visualization_msgs.msg import *
import rospy, cv2
from fast_possion import fast_poisson




class heightmap_recover:

    def __init__(self):
        self.kernal = self.make_kernal(21)
        self.kernal1 = self.make_kernal(10)
        self.kernal2 = self.make_kernal(18)
        self.kernal3 = self.make_kernal(26)
        self.kernal4 = self.make_kernal(2)
        self.kernal5 = self.make_kernal(8)
        self.kernal6 = self.make_kernal(3)
        self.index = 0
        self.M = np.load('M_GS2.npy')
        self.gradx_table = np.load('lookup_table_gx.npy')
        self.grady_table = np.load('lookup_table_gy.npy')
        self.lowbar = 10
        self.highbar = 70
        self.cols = 480
        self.rows = 640
        self.xv, self.yv = np.meshgrid(np.linspace(0, 434, 435, endpoint = True),np.linspace(0, 469, 470, endpoint = True))
        self.red_range = [25,255]
        self.green_range = [135,255]
        self.blue_range = [30,110]
        self.image_sub = rospy.Subscriber("/rpi/gelsight/raw_image2/compressed",CompressedImage,self.call_back)


    #image processing
    def rgb2gray(self,rgb):
        return np.dot(rgb[...,:3], [0.299, 0.587, 0.114])

    def make_kernal(self,n):
        kernal = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(n,n))
        return kernal 

    def calibration_v2(self,img):
        imgw = cv2.warpPerspective(img, self.M, (self.rows, self.cols))
        # im_temp = imgw*np.dstack((self.ROImask,self.ROImask,self.ROImask))
        # imgwc = im_temp[10:,60:572]
        imgwc = imgw[10:,95:-110,:]
        im_cal = imgwc/self.img_blur*100
        return im_cal,imgwc

    def creat_mask(self,im_cal):
        img_gray = self.rgb2gray(im_cal).astype(np.uint8)
        ret,thresh1 = cv2.threshold(img_gray,self.thre,255,cv2.THRESH_BINARY)
        final_image2 = cv2.erode(thresh1, self.kernal4, iterations=1)
        final_image = cv2.dilate(final_image2, self.kernal5, iterations=1)
        return final_image
    
    def color2position_array(self,colorimage):
        mask_red1 = (colorimage[:,:,2] >= self.red_range[1]).astype(np.uint8)*int((self.red_range[1] - self.red_range[0])/2)
        mask_red2 = ((colorimage[:,:,2] < self.red_range[1]).astype(np.uint8))*((colorimage[:,:,2] >= self.red_range[0]).astype(np.uint8))*((colorimage[:,:,2] - self.red_range[0])/2).astype(int)
        x = mask_red1 + mask_red2

        mask_green1 = (colorimage[:,:,1] >= self.green_range[1]).astype(np.uint8)*int((self.green_range[1] - self.green_range[0])/2)
        mask_green2 = ((colorimage[:,:,1] < self.green_range[1]).astype(np.uint8))*((colorimage[:,:,1] >= self.green_range[0]).astype(np.uint8))*(((colorimage[:,:,1] - self.green_range[0])/2).astype(int))
        y = mask_green1 + mask_green2
        
        mask_blue1 = (colorimage[:,:,0] >= self.blue_range[1]).astype(np.uint8)*int((self.blue_range[1] - self.blue_range[0])/2)
        mask_blue2 = ((colorimage[:,:,0] < self.blue_range[1]).astype(np.uint8))*((colorimage[:,:,0] >= self.blue_range[0]).astype(np.uint8))*(((colorimage[:,:,0] - self.blue_range[0])/2).astype(int))
        z = mask_blue1 + mask_blue2
        return x,y,z
                
    def contact_detection(self,im):
        im_sub = im/self.img_blur_ref.astype(np.float32)*70
        im_gray = self.rgb2gray(im_sub).astype(np.uint8)

#        mask_brightness = im_gray < 75
        # cv2.imshow('contact_image',mask_brightness.astype(np.uint8)*255)
        # cv2.waitKey(0)
        im_canny = cv2.Canny(im_gray,self.lowbar,self.highbar)
        im_canny = im_canny * self.amask
#        cv2.imshow('edge_image',im_canny) 
#        cv2.waitKey(100)
        img_d = cv2.dilate(im_canny, self.kernal1, iterations=1)
        img_e = cv2.erode(img_d, self.kernal1, iterations=1)
        img_ee = cv2.erode(img_e, self.kernal2, iterations=1)
        contact = cv2.dilate(img_ee, self.kernal3, iterations=1).astype(np.uint8)
        pad = 20
        contact[:pad,:] = 0
        contact[-pad:,:] = 0
        contact[:,:pad] = 0
        contact[:,-pad:] = 0
        return contact

        
    def call_back(self,data):
        t = time.time()
        np_arr = np.fromstring(data.data, np.uint8)
        image_np = cv2.imdecode(np_arr, cv2.IMREAD_COLOR)
        img = cv2.flip(image_np, 0) 
        
        if self.index < 10:
            imgw = cv2.warpPerspective(img, self.M, (self.rows,self.cols))
            imgwc = imgw[10:,95:-110,:]
            self.img_blur_ref = cv2.GaussianBlur(imgwc.astype(np.float32),(31,31),30)
            im_cal = imgwc/self.img_blur_ref*70
            im_gray = self.rgb2gray(im_cal).astype(np.uint8)
            self.amask = im_gray < 73
            self.index += 1
        else:
                    
            imgw = cv2.warpPerspective(img, self.M, (self.rows,self.cols))
            im_test = imgw[10:,95:-110,:]
            im_cal = im_test/self.img_blur_ref*100
#            cv2.imshow('test',im_test)
            img_gray = self.rgb2gray(im_cal).astype(np.uint8)
            ret,thresh1 = cv2.threshold(img_gray,75,255,cv2.THRESH_BINARY)
            final_image3 = cv2.dilate(thresh1,self.kernal4, iterations=1)
            final_image2 = cv2.erode(final_image3, self.kernal5, iterations=1)
            final_image = cv2.dilate(final_image2,self.kernal6, iterations=1)
            
            final_mask_color = np.concatenate((np.expand_dims(final_image,2),np.expand_dims(final_image,2),np.expand_dims(final_image,2)),axis = 2)
            img_blur = cv2.GaussianBlur(im_test,(31,31),30)
            im_final = (1-final_mask_color/255)*img_blur + im_test*(final_mask_color/255)
            img_lowpass = cv2.GaussianBlur(im_final,(9,13),2.5).astype(np.uint8)
            #cv2.imshow('low_pass_image',img_lowpass.astype(np.uint8))
            #cv2.waitKey(100)
            
            contact = self.contact_detection(im_test)
            contact_region = np.concatenate((np.expand_dims(contact,2),np.expand_dims(contact,2),np.expand_dims(contact,2)),axis = 2)/255
            img_lowpass = (img_lowpass*contact_region).astype(np.uint8)
            
            #im2show = (im_test*contact_region).astype(np.uint8)
            #plt.imshow(im2show)
            #im2show = (img_lowpass).astype(np.uint8)
            #cv2.imshow('low_pass',im2show)
            
            im_dx = np.zeros(img_lowpass[:,:,0].shape)
            im_dy = im_dx.copy()
            x_array, y_array, z_array = self.color2position_array(img_lowpass)
            im_dx = self.gradx_table[x_array,y_array,z_array]
            im_dy = self.grady_table[x_array,y_array,z_array]
            
            im_dx_smaller = cv2.pyrDown(im_dx)
            im_dy_smaller = cv2.pyrDown(im_dy)
            depth = fast_poisson(im_dy_smaller,im_dx_smaller)/im_dx_smaller.shape[0]/im_dx_smaller.shape[1]
#            depth2show = ((depth-np.amin(depth))/(np.amax(depth)+np.amin(depth))*100).astype(np.uint8)
            depth2show = ((depth+1)*150).astype(np.uint8)
#            depth = (depth>0)*depth
            depthC = cv2.applyColorMap(depth2show, cv2.COLORMAP_JET)
            cv2.imshow('rawimage',im_test.astype(np.uint8))
            cv2.imshow('heightmap',depthC.astype(np.uint8))
            key = cv2.waitKey(1)
            if key == '133':
                self.index = 1
                print 'initial frame reset'
            print 1/(time.time()-t)

            
def main():
    print "start"
    hr = heightmap_recover()
    rospy.init_node('heightmap', anonymous=True)
    rospy.spin()

if __name__ == "__main__": 
    main()
    
#%%
    
    
