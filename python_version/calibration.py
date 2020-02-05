import cv2
import numpy as np 
import matplotlib.pyplot as plt  







if __name__=="__main__":
    ref_img = cv2.imread('./new_data/ref.jpg')
    cv2.imshow('ref_image', ref_img)
    cv2.waitKey(0)
