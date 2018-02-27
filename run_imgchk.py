#!/usr/bin/python
# ******************************************************************************
# Name:    Image Subset Checker
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import sys
pyver_tup   = sys.version_info[0:2]
pyver_major = sys.version_info.major
pyver_minor = sys.version_info.minor

if pyver_tup >= (3, 5):
    pass
elif pyver_tup >= (2, 7) and pyver_tup <= (3,0):
    pass
else:
    raise Exception("Minimum Python ver 3.5 or 2.7 are required")

import os
# from PIL import Image
import cv2
import fx_img_process
import numpy as np
from scipy.misc import imread # uses PIL

# from PIL import Image

#>******************************************************************************
def main(imgpath_x, imgpath_y):
    if checkfile(imgpath_x):
        if checkfile(imgpath_y):
            fx_img_process.get_imgdata(imgpath_x, imgpath_y)
            # tempargs = base_procsort(imgpath_x, imgpath_y)
            # print(tempargs)
        else:
            sys.exit("Image file DNE")
    else:
        sys.exit("Image file DNE")

    print("Python %s.%s" % (pyver_major, pyver_minor))


#>******************************************************************************
def base_procsort(imgx, imgy):

    # ix = Image.open(imgx)
    # iy = Image.open(imgy)

    # ixformat = ix.format
    # ixsize = ix.size
    # ixtype = ix.mode
    # print("Size: %d %d Type: %s Source: %s" % (ixsize[0], ixsize[1], ixtype, ixformat))
    # iyformat = iy.format
    # iysize = iy.size
    # iytype = iy.mode

    # reimgx = cv2.imread(imgx)
    # reimgy = cv2.imread(imgy)
    reimgx = imread(imgx)
    reimgy = imread(imgy)

    hx, wx, chx = np.shape(reimgx)
    print("Height: %s" % str(hx))
    print("Width: %s" % str(wx))
    print("BPP : %s " % str(chx))

    hy, wy, chy = np.shape(reimgy)
    print("Height: %s" % str(hy))
    print("Width: %s" % str(wy))
    print("BPP : %s " % str(chy))

    if int(hx*wx) > int(hy*wy):
        return imgx, imgy
    elif int(hx*wx) < int(hy*wy):
        return imgy, imgx
    else:
        return 0
        '''How are the two images == in size?! '''


#>******************************************************************************
def checkfile(ifilepath):
    """ Checks if file exists. Follows symlinks """
    if os.path.isfile(ifilepath):
        return 1
    else:
        return 0



if __name__ == "__main__":
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        sys.exit("Two Terminal arguments are required")
