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

    if checkfile(imgpath_x):
        if checkfile(imgpath_y):
            areax, xrow, xcol = ret_imgarea(imgpath_x)
            areay, yrow, ycol = ret_imgarea(imgpath_y)
            imgnamex = os.path.basename(imgpath_x)
            imgnamey = os.path.basename(imgpath_y)
            print("Checking %s and %s" % (imgnamex, imgnamey))
            if areax>areay and xrow>=yrow and xcol>ycol or (areax>areay and xrow>yrow and xcol>=ycol):
                res_tup, final_rc, stat_str = fx_img_process.fft2_crosscorr(imgpath_x, imgpath_y)
                imgtype = res_tup[0]
                if stat_str=="True":
                    if imgtype=="Orig":
                        print("%s is a subimage of %s " % (imgnamey, imgnamex))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                    elif imgtype=="Flip":
                        print("%s is a flipped subimage of %s " % (imgnamey, imgnamex))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                elif stat_str=="False":
                    print(" %s is not a subimage of %s" % (imgnamey, imgnamex))
            elif areax<areay and xrow<=yrow and xcol<ycol or (areax<areay and xrow<yrow and xcol<=ycol):
                res_tup, final_rc, stat_str = fx_img_process.fft2_crosscorr(imgpath_y, imgpath_x)
                imgtype = res_tup[0]
                if stat_str=="True":
                    if imgtype=="Orig":
                        print("%s is a subimage of %s " % (imgnamex, imgnamey))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                    elif imgtype=="Flip":
                        print(" %s is a a flipped subimage of %s " % (imgnamex, imgnamey))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                elif stat_str=="False":
                    print(" %s is not a subimage of %s" % (imgnamex, imgnamey))
            else:
                sys.exit("Images uncomparable due to mismatched row/column dimensions") 
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
            imgnamex = os.path.basename(image_lt[x[0]])
            imgnamey = os.path.basename(image_lt[x[1]])
            if area_lt[x[0]] > area_lt[x[1]]:
                print("%s is larger." % imgnamex)
                # results_dict = fx_img_process.fft2_crosscorr(image_lt[x[0]], image_lt[x[1]])
                res_tup, final_rc, stat_str = fx_img_process.fft2_crosscorr(image_lt[x[0]], image_lt[x[1]])
                imgtype = res_tup[0]
                # print("imgtype: ", imgtype)
                if stat_str=="True":
                    if imgtype=="Orig":
                        print("%s is a subimage of %s " % (imgnamey, imgnamex))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                    elif imgtype=="Flip":
                        print("%s is a flipped subimage of %s " % (imgnamey, imgnamex))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                elif stat_str=="False":
                    print("%s is not a subimage of %s" % (imgnamey, imgnamex))
            elif area_lt[x[0]] < area_lt[x[1]]:
                print("%s is larger." % imgnamey)
                # results_dict = fx_img_process.fft2_crosscorr(image_lt[x[1]], image_lt[x[0]])
                res_tup, final_rc, stat_str = fx_img_process.fft2_crosscorr(image_lt[x[1]], image_lt[x[0]])
                imgtype = res_tup[0]
                # print("imgtype: ", imgtype)
                if stat_str=="True":
                    if imgtype=="Orig":
                        print("%s is a subimage of %s " % (imgnamex, imgnamey))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                    elif imgtype=="Flip":
                        print("%s is a a flipped subimage of %s " % (imgnamex, imgnamey))
                        print("Located at %i,%i  (row,col)\n" % (final_rc[0], final_rc[1]))
                elif stat_str=="False":
                    print("%s is not a subimage of %s" % (imgnamex, imgnamey))
            elif area_lt[x[0]] == area_lt[x[1]]:
                print("same size...")


#>******************************************************************************
def ret_imgarea(imgarg):
    """ Use scipy imread to get flattened/grayscale image data and calculate
    image area using row-col values from .shape """
    # img_dat = imread(imgpath, flatten=True)
    irow, icol = 0, 0
    if isinstance(imgarg, str):
        img_dat = ndimage.imread(imgarg, flatten=True)
        irow,icol = img_dat.shape
    else:
        irow, icol = imgarg.shape
    return irow*icol, irow, icol

#>******************************************************************************
def test_imgfeed(testpath=""):
    timgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/Testing/"
    timg1 = timgpath + 'jesusc8.jpg'
    timg2 = timgpath + 'jc8slice8.jpg'
    timg3 = timgpath + 'jc8slice8cut.jpg'
    timg4 = timgpath + 'jc8piece.jpg'
    tdiff1 = timgpath + 'wpKMy-minic.jpg'

    # return [timg1, timg2, timg3, timg4]
    return [[timg1, timg2], [timg1, timg3], [timg1, timg4], [timg1, tdiff1]]

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
    # print(lgimgpath)
    # lgimgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/LgTesting/"
    limg1 = "04041_mountrainier_2880x1800.jpg"
    limg1c = "04041_mountrainier_smcut.jpg"
    limg2 = "AJkBi5n.jpg"
    limg2c = "AJkBi5n_smcut.jpg"
    img_lt = [lgimgpath + xg for xg in [limg1, limg1c]]
    img_lt2 = [lgimgpath + xg for xg in [limg2, limg2c]]
    return [img_lt, img_lt2]

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
        try:
            test_main(sys.argv[1])
        except IndexError:
            timg_lt = test_imgfeed()
            if len(timg_lt) > 0:
                if isinstance(timg_lt[0], list):
                    for lt in timg_lt:
                        test_main(lt)
            else:
                sys.exit("Two Terminal arguments are required")
