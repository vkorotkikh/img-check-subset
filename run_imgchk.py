#!/usr/bin/python
# ******************************************************************************
# Name:    Image Subset Checker
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import os, sys
pyver_tup   = sys.version_info[0:2]
pyver_major = sys.version_info.major
pyver_minor = sys.version_info.minor

if pyver_tup >= (3, 5):
    pass
# elif pyver_tup >= (2, 7) and pyver_tup <= (3,0):
#     pass
else:
    # raise Exception("Minimum Python ver 3.5 or 2.7 are required")
    raise Exception("Minimum Python version 3.5 required")
    sys.exit()

import os
# from PIL import Image
import itertools as itr
import numpy as np
from scipy import fftpack, ndimage
from scipy.misc import imread # uses PIL
import fx_img_process

# from PIL import Image

#>******************************************************************************
def main(imgpath_x, imgpath_y):
    """ Check if files exist and then pass images ordered by size """
    rowinds, colinds = 0, 0
    if checkfile(imgpath_x):
        if checkfile(imgpath_y):
            areax, areay = ret_imgarea(imgpath_x), ret_imgarea(imgpath_y)
            if areax > areay:
                results_dict = fx_img_process.fft2_crosscorr(imgpath_x, imgpath_y)
                # fx_img_process.fft2_crosscorr(imgpath_x, imgpath_y)
            elif areax < areay:
                results_dict = fx_img_process.fft2_crosscorr(imgpath_y, imgpath_x)
                # fx_img_process.fft2_crosscorr(imgpath_y, imgpath_x)
            else:
                pass
        else:
            sys.exit("Image file DNE")
    else:
        sys.exit("Image file DNE")

#>******************************************************************************
def do_imgtrueflat(ipath):
    """ Using scipy.ndimage parse in image, flatten the data to grayscale and
    return grayscale 2d data """
    idata = ndimage.imread(ipath, flatten=True)
    return idata

#>******************************************************************************
def do_imgdata(ipath):
    """ Using scipy.ndimage to parse in image as RGB and return this data """
    return ndimage.imread(ipath, mode='RGB')

#>******************************************************************************
def do_imgzncc(ipath):
    """ normalize the data before doing cross-correlation calculations """
    imgdata = do_imgtrueflat(ipath)
    zndata = (imgdata - imgdata.mean())/ imgdata.std()
    return zndata

#>******************************************************************************
def do_imgfft2(ipath):
    imgdata = do_imgtrueflat(ipath)
    fft2dat = fftpack.fft2(imgdata)
    return fft2dat

#>******************************************************************************
def test_main(image_lt):
    """ main used for testing purposes. Uses hardcoded image paths and names """
    # narg = [do_imgzncc(ix) for ix in image_lt]
    ifdat_lt = [do_imgzncc(ix) for ix in image_lt]
    area_lt = [ret_imgarea(ix) for ix in ifdat_lt]
    acombs = list(itr.combinations(list(range(0,len(image_lt))),2))

    rowinds, colinds = 0, 0
    for x in acombs:
        if x[0] == 0 or x[0] == 1:
            if area_lt[x[0]] > area_lt[x[1]]:
                print("Img 1 is bigger")
                results_dict = fx_img_process.fft2_crosscorr(image_lt[x[0]], image_lt[x[1]])
                # rowinds, colinds = fx_img_process.fft2_crosscorr(image_lt[x[0]], image_lt[x[1]])
            elif area_lt[x[0]] < area_lt[x[1]]:
                print("Img 2 is bigger")
                results_dict = fx_img_process.fft2_crosscorr(image_lt[x[1]], image_lt[x[0]])
                # rowinds, colinds = fx_img_process.fft2_crosscorr(image_lt[x[1]], image_lt[x[0]])
            elif area_lt[x[0]] == area_lt[x[1]]:
                print("same size...")


#>******************************************************************************
def ret_imgarea(imgarg):
    """ Use scipy imread to get flattened/grayscale image data and calculate
    image area using row-col values from .shape """
    # img_dat = imread(imgpath, flatten=True)
    irow, icol = imgarg.shape
    return irow*icol

#>******************************************************************************
def test_imgfeed(testpath=""):
    timgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/Testing/"
    timg1 = timgpath + 'jesusc8.jpg'
    timg2 = timgpath + 'jc8slice8.jpg'
    timg3 = timgpath + 'jc8slice8cut.jpg'
    timg4 = timgpath + 'jc8piece.jpg'
    tdiff1 = timgpath + 'wpKMy-minic.jpg'

    # return [timg1, timg2, timg3, timg4]
    return [timg1, timg3]

#>******************************************************************************
def test_bigimg(testpath=""):
    """ Return list of all large images and their subimages/slices
    testpath - Maybe(?) implemented later
    """
    exefilepath = os.path.realpath(__file__)
    exedirpath = os.path.dirname(exefilepath)
    lgimgpath = ""

    if os.path.isdir((exedirpath + "/LgTesting")):
        lgimgpath =  exedirpath + "/LgTesting/"
    print(lgimgpath)
    # lgimgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/LgTesting/"
    limg1 = "04041_mountrainier_2880x1800.jpg"
    limg1c = "04041_mountrainier_smcut.jpg"
    limg2 = "AJkBi5n.jpg"
    limg2c = "AJkBi5n_smcut.jpg"
    img_lt = [lgimgpath + xg for xg in [limg1, limg1c]]
    img_lt2 = [lgimgpath + xg for xg in [limg2, limg2c]]
    return img_lt2

#>******************************************************************************
def checkfile(ifilepath):
    """ Checks if file exists at given path=ifilepath. Follows symlinks """
    if os.path.isfile(ifilepath):
        return 1
    else:
        return 0


#>******************************************************************************
if __name__ == "__main__":
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        timg_lt = test_bigimg()
        if len(timg_lt) > 0:
            test_main(timg_lt)
        else:
            sys.exit("Two Terminal arguments are required")
