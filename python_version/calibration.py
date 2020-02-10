import cv2
import numpy as np 
import matplotlib.pyplot as plt  
import glob
from fast_poisson import fast_poisson
#from fast_poisson_shawn import poisson_reconstruct
from mpl_toolkits.mplot3d import Axes3D
from scipy import signal
from skimage.restoration import inpaint

class image_processor:
    def __init__(self):
        pass
    
    def crop_image(self,img, pad):
        return img[pad:-pad,pad:-pad]

class calibration:
    def __init__(self):
        self.BallRad=6.35/2 #mm
        self.Pixmm = 0.0806
        self.ratio = 1/2.
        self.red_range = [-90, 90]
        self.green_range = [-90, 90] #[-60, 50]
        self.blue_range = [-90, 90] # [-80, 60]
        self.red_bin = int((self.red_range[1] - self.red_range[0])*self.ratio)
        self.green_bin = int((self.green_range[1] - self.green_range[0])*self.ratio)
        self.blue_bin = int((self.blue_range[1] - self.blue_range[0])*self.ratio)
        self.zeropoint = -90;
        self.lookscale = 180;
        self.bin_num = 90
    
    def mask_marker(self, raw_image):
        m, n = raw_image.shape[1], raw_image.shape[0]
        raw_image = cv2.pyrDown(raw_image).astype(np.float32)
        blur = cv2.GaussianBlur(raw_image, (25, 25), 0)
        blur2 = cv2.GaussianBlur(raw_image, (5, 5), 0)
        diff = blur - blur2
        diff *= 16.0
        # cv2.imshow('blur2', blur.astype(np.uint8))
        # cv2.waitKey(1)

        diff[diff < 0.] = 0.
        diff[diff > 255.] = 255.

        # diff = cv2.GaussianBlur(diff, (5, 5), 0)
        # cv2.imshow('diff', diff.astype(np.uint8))
        # cv2.waitKey(1)
        mask = (diff[:, :, 0] > 25) & (diff[:, :, 2] > 25) & (diff[:, :, 1] >
                                                              120)
        # cv2.imshow('mask', mask.astype(np.uint8) * 255)
        # cv2.waitKey(1)
        mask = cv2.resize(mask.astype(np.uint8), (m, n))
#        mask = mask * self.dmask
#        mask = cv2.dilate(mask, self.kernal4, iterations=1)

        # mask = cv2.erode(mask, self.kernal4, iterations=1)
        return (1 - mask) * 255
    
    def find_dots(self, binary_image):
        # down_image = cv2.resize(binary_image, None, fx=2, fy=2)
        params = cv2.SimpleBlobDetector_Params()
        # Change thresholds
        params.minThreshold = 1
        params.maxThreshold = 12
        params.minDistBetweenBlobs = 9
        params.filterByArea = True
        params.minArea = 9
        params.filterByCircularity = False
        params.filterByConvexity = False
        params.filterByInertia = False
        params.minInertiaRatio = 0.5
        detector = cv2.SimpleBlobDetector_create(params)
        keypoints = detector.detect(binary_image.astype(np.uint8))
        # im_to_show = (np.stack((binary_image,)*3, axis=-1)-100)
        # for i in range(len(keypoints)):
        #     cv2.circle(im_to_show, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 5, (0, 100, 100), -1)

        # cv2.imshow('final_image1',im_to_show)
        # cv2.waitKey(1)
        return keypoints
    
    def make_mask(self,img, keypoints):
        img = np.zeros_like(img[:,:,0])
        for i in range(len(keypoints)):
            cv2.circle(img, (int(keypoints[i].pt[0]), int(keypoints[i].pt[1])), 6, (1), -1)

        return img
    
    def contact_detection(self,raw_image, ref, marker_mask):
        blur = cv2.GaussianBlur(ref.astype(np.float32), (25, 25), 0)
        diff_img = np.max(np.abs(raw_image.astype(np.float32) - blur),axis = 2)
        contact_mask = (diff_img> 30).astype(np.uint8)*(1-marker_mask)
        contours,_ = cv2.findContours(contact_mask,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        areas = [cv2.contourArea(c) for c in contours]
        sorted_areas = np.sort(areas)
        cnt=contours[areas.index(sorted_areas[-1])] #the biggest contour
        (x,y),radius = cv2.minEnclosingCircle(cnt)
        center = (int(x),int(y))
        radius = int(radius)
        
        key = -1
        while key != 27:
            center = (int(x),int(y))
            radius = int(radius)
            im2show = cv2.circle(np.array(raw_image),center,radius,(0,40,0),2)
            cv2.imshow('contact', im2show.astype(np.uint8))
            key = cv2.waitKey(0)
            if key == 119:
                y -= 1
            elif key == 115:
                y += 1
            elif key == 97:
                x -= 1
            elif key == 100:
                x += 1
            elif key == 109:
                radius += 1
            elif key == 110:
                radius -= 1

        contact_mask = np.zeros_like(contact_mask)
        cv2.circle(contact_mask,center,radius,(1),-1)
        contact_mask = contact_mask * (1-marker_mask)
#        cv2.imshow('contact_mask',contact_mask*255)
#        cv2.waitKey(0)
        return contact_mask, center, radius

    def get_gradient(self, img, ref, center, radius_p, valid_mask, table, table_account):
        ball_radius_p = self.BallRad / self.Pixmm
        blur = cv2.GaussianBlur(ref.astype(np.float32), (25, 25), 0)
        img_smooth = cv2.GaussianBlur(img.astype(np.float32), (3, 3), 0)
        diff = img_smooth - blur 
#        diff_valid = np.abs(diff * np.dstack((valid_mask,valid_mask,valid_mask)))
        pixels_valid = diff[valid_mask>0]
        pixels_valid[:,0] = np.clip((pixels_valid[:,0] - self.blue_range[0])*self.ratio, 0, self.blue_bin-1)
        pixels_valid[:,1] = np.clip((pixels_valid[:,1] - self.green_range[0])*self.ratio, 0, self.green_bin-1)
        pixels_valid[:,2] = np.clip((pixels_valid[:,2] - self.red_range[0])*self.ratio, 0, self.red_bin-1)
        pixels_valid = pixels_valid.astype(int)
        
        
#        range_blue = [max(np.mean(pixels_valid[:,0])-2*np.std(pixels_valid[:,0]),np.min(pixels_valid[:,0])), \
#                     min(np.mean(pixels_valid[:,0])+2.5*np.std(pixels_valid[:,0]), np.max(pixels_valid[:,0]))] 
#        range_green = [max(np.mean(pixels_valid[:,1])-2*np.std(pixels_valid[:,1]),np.min(pixels_valid[:,1])), \
#                     min(np.mean(pixels_valid[:,1])+2.5*np.std(pixels_valid[:,1]), np.max(pixels_valid[:,1]))] 
#        range_red = [max(np.mean(pixels_valid[:,2])-2*np.std(pixels_valid[:,2]),np.min(pixels_valid[:,2])), \
#                     min(np.mean(pixels_valid[:,2])+2.5*np.std(pixels_valid[:,2]), np.max(pixels_valid[:,2]))] 
        
#        print('blue', range_blue, 'green', range_green, 'red',  range_red)
        
#        print(np.min(pixels_valid[:,0]), np.max(pixels_valid[:,0]))
#        print(np.min(pixels_valid[:,1]), np.max(pixels_valid[:,1]))
#        print(np.min(pixels_valid[:,2]), np.max(pixels_valid[:,2]))
#        print(pixels_valid.shape)
#        plt.figure(0)
#        plt.hist(pixels_valid[:,0], bins = 256)
#        plt.figure(1)
#        plt.hist(pixels_valid[:,1], bins = 256)
#        plt.figure(2)
#        plt.hist(pixels_valid[:,2], bins = 256)
#        plt.show()
        x = np.linspace(0, img.shape[0],img.shape[0])
        y = np.linspace(0, img.shape[1],img.shape[1])
        xv, yv = np.meshgrid(y, x)
#        print('img shape', img.shape, xv.shape, yv.shape)
        xv = xv - center[0]
        yv = yv - center[1]
        rv = np.sqrt(xv**2 + yv**2)
        radius_p = min(radius_p, ball_radius_p-1) 
        mask = (rv < radius_p)
        mask_small = (rv < radius_p-1)
#        gradmag=np.arcsin(rv*mask/ball_radius_p)*mask;
#        graddir=np.arctan2(-yv*mask, -xv*mask)*mask;
#        gradx_img=gradmag*np.cos(graddir);
#        grady_img=gradmag*np.sin(graddir);
#        depth = fast_poisson(gradx_img, grady_img)
        temp = ((xv*mask)**2 + (yv*mask)**2)*self.Pixmm**2
        height_map = (np.sqrt(self.BallRad**2-temp)*mask - np.sqrt(self.BallRad**2-(radius_p*self.Pixmm)**2))*mask
        height_map[np.isnan(height_map)] = 0
#        depth = poisson_reconstruct(grady_img, gradx_img, np.zeros(grady_img.shape))
        gx_num = signal.convolve2d(height_map, np.array([[0,0,0],[0.5,0,-0.5],[0,0,0]]), boundary='symm', mode='same')*mask_small
        gy_num = signal.convolve2d(height_map, np.array([[0,0,0],[0.5,0,-0.5],[0,0,0]]).T, boundary='symm', mode='same')*mask_small
#        depth_num = fast_poisson(gx_num, gy_num)
#        plt.imshow(height_map)
#        plt.show()
        gradxseq = gx_num[valid_mask>0]
        gradyseq = gy_num[valid_mask>0]
        
        for i in range(gradxseq.shape[0]):
            b, g, r = pixels_valid[i,0], pixels_valid[i,1], pixels_valid[i,2]
#            print(r,g,b)
            if table_account[b,g,r] < 1.: 
                table[b,g,r,0] = gradxseq[i]
                table[b,g,r,1] = gradyseq[i]
                table_account[b,g,r] += 1
            else:
#                print(table[b,g,r,0], gradxseq[i], table[b,g,r,1], gradyseq[i])
                table[b,g,r,0] = (table[b,g,r,0]*table_account[b,g,r] + gradxseq[i])/(table_account[b,g,r]+1)
                table[b,g,r,1] = (table[b,g,r,1]*table_account[b,g,r] + gradyseq[i])/(table_account[b,g,r]+1)
                table_account[b,g,r] += 1
        return table, table_account
        
    

#        plt.figure(0) 
#        plt.imshow(gx_num)
#        plt.figure(1)
#        plt.imshow(gy_num)
#        plt.show()
      
    def get_gradient_v2(self, img, ref, center, radius_p, valid_mask, table, table_account):
        ball_radius_p = self.BallRad / self.Pixmm
        blur = cv2.GaussianBlur(ref.astype(np.float32), (49, 49), 49)
        blur_inverse = 1+((np.mean(blur)/blur)-1)*2;
        img_smooth = cv2.GaussianBlur(img.astype(np.float32), (3, 3), 0)
        diff_temp1 = img_smooth - blur 
        diff_temp2 = diff_temp1 * blur_inverse
        diff_temp3 = np.clip((diff_temp2-self.zeropoint)/self.lookscale,0,0.999)
        diff = (diff_temp3*self.bin_num).astype(int)
#        diff_valid = np.abs(diff * np.dstack((valid_mask,valid_mask,valid_mask)))
        pixels_valid = diff[valid_mask>0]
#        pixels_valid[:,0] = np.clip((pixels_valid[:,0] - self.blue_range[0])*self.ratio, 0, self.blue_bin-1)
#        pixels_valid[:,1] = np.clip((pixels_valid[:,1] - self.green_range[0])*self.ratio, 0, self.green_bin-1)
#        pixels_valid[:,2] = np.clip((pixels_valid[:,2] - self.red_range[0])*self.ratio, 0, self.red_bin-1)
#        pixels_valid = pixels_valid.astype(int)
        
        
#        range_blue = [max(np.mean(pixels_valid[:,0])-2*np.std(pixels_valid[:,0]),np.min(pixels_valid[:,0])), \
#                     min(np.mean(pixels_valid[:,0])+2.5*np.std(pixels_valid[:,0]), np.max(pixels_valid[:,0]))] 
#        range_green = [max(np.mean(pixels_valid[:,1])-2*np.std(pixels_valid[:,1]),np.min(pixels_valid[:,1])), \
#                     min(np.mean(pixels_valid[:,1])+2.5*np.std(pixels_valid[:,1]), np.max(pixels_valid[:,1]))] 
#        range_red = [max(np.mean(pixels_valid[:,2])-2*np.std(pixels_valid[:,2]),np.min(pixels_valid[:,2])), \
#                     min(np.mean(pixels_valid[:,2])+2.5*np.std(pixels_valid[:,2]), np.max(pixels_valid[:,2]))] 
        
#        print('blue', range_blue, 'green', range_green, 'red',  range_red)
        
#        print(np.min(pixels_valid[:,0]), np.max(pixels_valid[:,0]))
#        print(np.min(pixels_valid[:,1]), np.max(pixels_valid[:,1]))
#        print(np.min(pixels_valid[:,2]), np.max(pixels_valid[:,2]))
#        print(pixels_valid.shape)
#        plt.figure(0)
#        plt.hist(pixels_valid[:,0], bins = 256)
#        plt.figure(1)
#        plt.hist(pixels_valid[:,1], bins = 256)
#        plt.figure(2)
#        plt.hist(pixels_valid[:,2], bins = 256)
#        plt.show()
        x = np.linspace(0, img.shape[0],img.shape[0])
        y = np.linspace(0, img.shape[1],img.shape[1])
        xv, yv = np.meshgrid(y, x)
#        print('img shape', img.shape, xv.shape, yv.shape)
        xv = xv - center[0]
        yv = yv - center[1]
        rv = np.sqrt(xv**2 + yv**2)
        radius_p = min(radius_p, ball_radius_p-1) 
        mask = (rv < radius_p)
        mask_small = (rv < radius_p-1)
#        gradmag=np.arcsin(rv*mask/ball_radius_p)*mask;
#        graddir=np.arctan2(-yv*mask, -xv*mask)*mask;
#        gradx_img=gradmag*np.cos(graddir);
#        grady_img=gradmag*np.sin(graddir);
#        depth = fast_poisson(gradx_img, grady_img)
        temp = ((xv*mask)**2 + (yv*mask)**2)*self.Pixmm**2
        height_map = (np.sqrt(self.BallRad**2-temp)*mask - np.sqrt(self.BallRad**2-(radius_p*self.Pixmm)**2))*mask
        height_map[np.isnan(height_map)] = 0
#        depth = poisson_reconstruct(grady_img, gradx_img, np.zeros(grady_img.shape))
        gx_num = signal.convolve2d(height_map, np.array([[0,0,0],[0.5,0,-0.5],[0,0,0]]), boundary='symm', mode='same')*mask_small
        gy_num = signal.convolve2d(height_map, np.array([[0,0,0],[0.5,0,-0.5],[0,0,0]]).T, boundary='symm', mode='same')*mask_small
#        depth_num = fast_poisson(gx_num, gy_num)
#        plt.imshow(height_map)
#        plt.show()
        gradxseq = gx_num[valid_mask>0]
        gradyseq = gy_num[valid_mask>0]
        
        for i in range(gradxseq.shape[0]):
            b, g, r = pixels_valid[i,0], pixels_valid[i,1], pixels_valid[i,2]
#            print(r,g,b)
            if table_account[b,g,r] < 1.: 
                table[b,g,r,0] = gradxseq[i]
                table[b,g,r,1] = gradyseq[i]
                table_account[b,g,r] += 1
            else:
#                print(table[b,g,r,0], gradxseq[i], table[b,g,r,1], gradyseq[i])
                table[b,g,r,0] = (table[b,g,r,0]*table_account[b,g,r] + gradxseq[i])/(table_account[b,g,r]+1)
                table[b,g,r,1] = (table[b,g,r,1]*table_account[b,g,r] + gradyseq[i])/(table_account[b,g,r]+1)
                table_account[b,g,r] += 1
        return table, table_account
        
    

#        plt.figure(0) 
#        plt.imshow(gx_num)
#        plt.figure(1)
#        plt.imshow(gy_num)
#        plt.show()   
        

if __name__=="__main__":
    cali = calibration()
    imp = image_processor()
    pad = 20
    ref_img = cv2.imread('./test_data/ref.jpg')
    ref_img = imp.crop_image(ref_img, pad)
    table = np.zeros((cali.blue_bin, cali.green_bin, cali.red_bin, 2))
    table_account = np.zeros((cali.blue_bin, cali.green_bin, cali.red_bin))
#    cv2.imshow('ref_image', ref_img)
#    cv2.waitKey(0)
    has_marke = True
    img_list = glob.glob("test_data/sample*.jpg")
    
    for name in img_list:
#        print(name)
        img = cv2.imread(name)
        img = imp.crop_image(img, pad)
        if has_marke: 
            marker = cali.mask_marker(img)
            keypoints = cali.find_dots(marker)
            marker_mask = cali.make_mask(img, keypoints)
#            im2show = img.copy()
#            im2show[:,:,1] += marker_mask*40
#            cv2.imshow('marker_mask', im2show.astype(np.uint8)) 
#            cv2.waitKey(0) 
        else:
            marker_mask = np.zeros_like(img)
        valid_mask, center, radius_p  = cali.contact_detection(img, ref_img, marker_mask)
        table, table_account = cali.get_gradient_v2(img, ref_img, center, radius_p, valid_mask, table, table_account)
    np.save('table_2.npy', table)
    np.save('count_map_2.npy', table_account)

#%%
#def make_kernal(n,k_type):
#    if k_type == 'circle':
#        kernal = cv2.getStructuringElement(cv2.MORPH_ELLIPSE,(n,n))
#    else:
#        kernal = cv2.getStructuringElement(cv2.MORPH_RECT,(n,n))
#    return kernal 
#
#table = np.load('table.npy')
#mask = (np.abs(table[:,:,:,0])>0).astype(np.uint8)
#mask_orginal = mask.copy()
#kernal = make_kernal(9,'circle')
#for i in range(90):
#    mask[:,:,i] = cv2.dilate(mask[:,:,i], kernal, iterations=1)
#    mask[:,:,i] = cv2.erode(mask[:,:,i], kernal, iterations=1)
#
#mask = mask - mask_orginal
#%%




#table_smooth = table.copy()
#
##for j in range(2):
##    for i in range(90):
#
#for i in range(90):
#    table_smooth[:,:,i,0] = inpaint.inpaint_biharmonic(table_smooth[:,:,i,0], \
#                mask[:,:,i],multichannel=False)
#    table_smooth[:,:,i,1] = inpaint.inpaint_biharmonic(table_smooth[:,:,i,1], \
#                mask[:,:,i],multichannel=False)
#    print(i)

   
#%%
#from scipy.interpolate import griddata
#table = np.load('table.npy')
#mask = (np.abs(table[:,:,:,0])>0).astype(np.uint8)
#grid_x, grid_y = np.mgrid[0:90, 0:90]
#table_smooth = table.copy()
#for i in range(90):
#    x_coor = grid_x[mask[:,:,i]>0]
#    y_coor = grid_y[mask[:,:,i]>0]
#    coor_data = np.vstack((x_coor, y_coor)).T
#    datax_temp = table_smooth[:,:,i,0]
#    datay_temp = table_smooth[:,:,i,1]
#    data_x = datax_temp[mask[:,:,i]>0]
#    data_y = datay_temp[mask[:,:,i]>0]
#    
##    print(coor_data.shape)
#    if data_x.shape[0] >0:
#        table_smooth[:,:,i,0] = griddata(coor_data, data_x, (grid_x, grid_y), method='linear')
#        table_smooth[:,:,i,1] = griddata(coor_data, data_y, (grid_x, grid_y), method='linear')
#    print(i/90.)
#    
#    
#
#table_smooth[np.isnan(table_smooth)] = 0.
#%%
#num = 30
#plt.figure(0)
#plt.imshow(table[:,:,num,0])
#plt.figure(1)
#plt.imshow(table_smooth[:,:,num,0])
#plt.figure(2)
#plt.imshow(mask[:,:,num])
#plt.show()










