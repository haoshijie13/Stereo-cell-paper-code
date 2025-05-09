import numpy as np
import skimage
import math
import pandas as pd
import sys
import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = pow(2,50).__str__()
import cv2 as cv
img = cv.imread(sys.argv[1], cv.IMREAD_GRAYSCALE)
assert img is not None, "file could not be read, check with os.path.exists()"
ret,thresh = cv.threshold(img,127,255,0)
contours,hierarchy = cv.findContours(thresh, 1, 2)


# Create an empty image for contours
img_contours = np.zeros(img.shape,np.uint32)

iter=0
stat=pd.DataFrame([],columns=["area","perimeter","Cx","Cy","AR","solidity","roundness","circularity","cell"])

for i in contours:
    area = cv.contourArea(i)
    perimeter = cv.arcLength(i,True)
    x,y,w,h = cv.boundingRect(i)
    ar = float(w)/h
    hull = cv.convexHull(i)
    hull_area = cv.contourArea(hull)
    if(hull_area==0):
        continue
    solidity = float(area)/hull_area
    roundness=4*area/(math.pi*w**2)
    circularity=4*math.pi*area/perimeter**2
    M = cv.moments(i)
    cx = round(M['m10']/M['m00'],ndigits=3)
    cy = round(M['m01']/M['m00'],ndigits=3)
    iter=iter+1
    #cv.fillPoly(img_contours, pts = [i], color=iter)
    rr, cc = skimage.draw.polygon(tuple(i[:,:,1].reshape(i.shape[0])),tuple(i[:,:,0].reshape(i.shape[0])))
    img_contours[rr,cc]=iter
    stat.loc[len(stat)] = [area,perimeter,cx,cy,ar,solidity,roundness,circularity,iter]

df = pd.DataFrame(np.vstack((np.where(img_contours>0)[0],np.where(img_contours>0)[1],img_contours[np.where(img_contours>0)])).T, columns = ['y','x','cell'])    
skimage.io.imsave('contours.png',img_contours)
df.to_csv(sys.argv[2]+".cell_mask.csv")
stat.to_csv(sys.argv[2]+".measurements.csv")
