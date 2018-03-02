#!/usr/bin/python
# ******************************************************************************
# Name:    Image Subset Checker
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import sys
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import itertools as itr
from pprint import pprint
from scipy import fftpack, ndimage, average
from scipy.misc import imread
from scipy.signal.signaltools import correlate2d as cor2d

#>******************************************************************************
def do_grayscale(imgarr):
    """ Simplistic grayscale operation on a imgdata array """
    if len(imgarr.shape) == 3:
        return sp.average(imgarr, -1)
    else:
        return imgarr

#>******************************************************************************
def do_imzncc(ipath):
    """ imread image data and perform normalization """
    imgdata = ndimage.imread(ipath, flatten=True)
    zndata = (imgdata - imgdata.mean())/ imgdata.std()
    return zndata

#>******************************************************************************
def do_imgfft2(ipath):
    imgdata = ndimage.imread(ipath, flatten=True)
    fft2dat = fftpack.fft2(imgdata)
    return fft2dat

#>******************************************************************************
def do_znccfft2(ipath):
    imgdata = ndimage.imread(ipath, flatten=True)
    zndata = (imgdata - imgdata.mean())/ imgdata.std()
    fft2dat = fftpack.fft2(zndata)
    return fft2dat

#>******************************************************************************
def doplt(ifftdata):
#     plt.imshow(20*np.log10(abs(ifftdata)))
    from matplotlib.colors import LogNorm
    plt.figure(figsize=(8, 6), dpi=80)
    plt.imshow(np.abs(ifftdata), norm=LogNorm(vmin=2))
#     plt.colorbar()
    plt.show()

#>******************************************************************************
def fft2_crosscorr(imgx_dat, imgy_dat):
    """ Assume imgx_dat is always going to be the bigger image, area wise"""
    imgx_dim = imgx_dat.shape
    print(imgx_dim)
    imgy_dim = imgy_dat.shape
    ixd, iyd = imgx_dim, imgy_dim
    xht = ixd[0] # heigth
    xwd = ixd[1] # width
    yht = iyd[0] # h
    ywd = iyd[1] # w

    divnum = 0
    mulfact = 4
    upperb, lowerb = 0, 0

    locmax_lt = []
    locnmax_lt = []
    if xht > yht and xwd == ywd: # loop over rows
        loc_lt = oned_slide_xcor(imgx_dat, imgy_dat, xht, yht, xwd, 'rows')
        for xloc in loc_lt:
            print(xloc)
#         divnum = xht//yht
#         stepd = mulfact*int(xht/divnum)
#         for i in range(0, divnum):
#             if (i*stepd + yht) < xht and i*(2*stepd) < xht:
#                 upperb = i*stepd + yht
#                 lowerb = i*stepd
#             else:
#                 upperb = xht
#                 lowerb = xht - yht
    elif xwd > ywd and xht == yht: # loop over columns
        loc_lt = oned_slide_xcor(imgx_dat, imgy_dat, xwd, ywd, xht, 'cols')
        for xloc in loc_lt:
            print(xloc)
    elif xht > yht and xwd > ywd:
        mulfact = 3
        divnumh = mulfact*(xht//yht) # Make stepsize 1/2 subimage size in that vec
        divnumw = mulfact*(xwd//ywd)
        # stepht = mulfact*int(xht/divnumh)
        stepht = int(xht/divnumh)
        stepwd = int(xwd/divnumw)

        upperht, lowerht = 0, 0
        upperwd, lowerwd = 0, 0

        subftconj = fftpack.fft2(imgy_dat).conj()
        print("Ht",xht,yht," Wd", xwd, ywd)
        print("Stepht: ", stepht, "Stephwd: ", stepwd)
        """ Build the index - pixel loc matrix """
        rcindices_lt = []
        for ih in range(0, divnumh):
            ''' Do stepsize everywhere, so stepht instead yht and ... '''
            tm_rcinds = []
            temp_lt = []
            tmin_lt = []
            if (ih*stepht + yht) < xht:
                upperht = ih*stepht + yht
                lowerht = ih*stepht
            else:
                upperht = xht
                lowerht = xht - yht
            rcindices_lt.append((lowerht, upperht))
            for iw in range(0, divnumw):
                if (iw*stepwd + ywd) < xwd:
#                 if (iw*stepwd + ywd) < (xwd-stepwd):
                    upperwd = iw*stepwd + ywd
                    lowerwd = iw*stepwd
                else:
                    upperwd = xwd
                    lowerwd = xwd - ywd
                tm_rcinds.append(((lowerht, upperht), (lowerwd, upperwd)))
                temp_dat = imgx_dat[lowerht:upperht, lowerwd:upperwd]
#                 print(temp_dat.shape, type(temp_dat))
#                 temp_dat = temp_dat[:, lowerwd:upperwd]
                tempd_min = (temp_dat - temp_dat.mean())/ temp_dat.std()
                locmax = twod_slide_xcorr(temp_dat, imgy_dat)
                locmaxnorm = twod_slide_xcorrfast(tempd_min, subftconj)
                rowstr = str(lowerht) + ":" + str(upperht)
                colstr = str(lowerwd) + ":" + str(upperwd)
                print("Indices", rowstr, colstr, " LocMax: ", locmax)
                temp_lt.append(int(locmax))
                tmin_lt.append(int(locmaxnorm))
            rcindices_lt.append(tm_rcinds)
            locmax_lt.append(temp_lt)
            # locnmax_lt.append(tmin_lt)
    locmax_arr = np.asarray(locmax_lt)
    print(locmax_arr.shape)
    print(locmax_arr)
    print(np.std(locmax_arr))
    locmax_minarr = (locmax_arr - int(locmax_arr.mean())) / (np.std(locmax_arr))
    locmax_minarr = locmax_minarr.astype(int)
    print("")
    print(np.mean(locmax_minarr), np.average(locmax_minarr))
    for xarr in locmax_minarr:
        print(xarr, "", np.average(xarr))
    for xlt in locmax_lt:
        print("Index", np.argmax(xlt), np.amax(xlt))
        print(xlt)
    print("")
    for xrc in rcindices_lt:
        print(xrc)
    locmax_slicing(locmax_minarr, rcindices_lt)
    # for ilt in locnmax_lt:
    #     print("Index", np.argmax(ilt), np.amax(ilt))
    #     print(ilt)

#>******************************************************************************
def locmax_slicing(lmax_array, rowcol_inds):
    """ Using the result values to cut up the original image into individual
    elements for more indepth subimage localization """
    print(type(lmax_array), "", lmax_array.shape)
    rows, cols = lmax_array.shape
    tmp_lmax = lmax_array.copy()

    newslice_lt = []
    tmp_hold = []

    for rn in range(0, rows):

        temp_s = []
        for cn in range(0, cols):
            if rn==0:
                if tmp_lmax[rn, cn] > 0:
                    temp_s.append(rowcol_inds[0][rn], rowcol_inds[1][cn])
                    if (cn+1) < cols and tmp_lmax[rn, cn+1] > 0:
                        pass
                    else:
                        tmp_hold.append(temp_s)
                        temp_s = []
                else:
                    pass
            elif rn>0:
                if tmp_lmax[rn, cn] > 0:
                    if tmp_lmax[rn-1, cn] > 0:
                        pass
                else:
                    pass


#>******************************************************************************
def twod_slide_xcorr(islice, subimg):
#     print("slice:", islice.shape, "sub:", subimg.shape)
    img_product = fftpack.fft2(islice) * fftpack.fft2(subimg).conj()
    inv_prod = fftpack.ifft2(img_product)
    return np.amax(inv_prod.real)
#     doplt(islice)

#>******************************************************************************
def twod_slide_xcorrfast(islice, simgfcon):
    img_product = fftpack.fft2(islice) * simgfcon
    inv_prod = fftpack.ifft2(img_product)
    return np.amax(inv_prod.real)

#>******************************************************************************
def oned_slide_xcor(mimg_dat, simg_dat, majx, minx, equald, ldirec):
    """ In the case when the subimage has either width or height = to that of
    big image
    Perform sliding 2d cross-correlate calculation over the smaller dimension
    mimg_dat - Main image numpy data array - normalized
    simg_dat - Subimage numpy data array - normalized
    majx - size major dimension of  full image - int
    minx - size minor dimension of subimage thats < majd
    equald - size of dimension thats = between full image and subimage      """
    mulfact = 4
    upperb, lowerb = 0, 0
    divnum = mulfact*(majx//minx)
    stepd = int(majx/divnum)

    cxcorr_maxlt = []
    for i in range(0, divnum):
        if i*stepd + minx < majx:
            upperb = i*stepd + minx
            lowerb = i*stepd
        else:
            upperb = majx
            lowerb = majx - minx
        if ldirec == 'rows':
            islc = mimg_dat[lowerb:upperb]
            print("Height:", islc.shape[0], "Width:", islc.shape[1])
            img_product = fftpack.fft2(islc) * fftpack.fft2(simg_dat).conj()
            inv_prod = fftpack.ifft2(img_product)
            cxcorr_maxlt.append((i, np.amax(inv_prod.real)))
#             doplt(islc)
#             doplt(inv_prod)
        elif ldirec == 'cols':
            islc = mimg_dat[:][lowerb:upperb]
            print("Height:", islc.shape[0], "Width:", islc.shape[1])
            img_product = fftpack.fft2(islc) * fftpack.fft2(simg_dat).conj()
            inv_prod = fftpack.fft2(img_product)
            cxcorr_maxlt.append((i, np.amax(inv_prod.real)))
    return cxcorr_maxlt


#>******************************************************************************
    # divn = 4*(ixd[0]//iyd[0])
    # stepd = int(ixd[0]/divn)
    #
    # for i in range(0, divn):
    #     highs = 0
    #     if i*stepd + iyd[0] < ixd[0]:
    #         highs = i*stepd + iyd[0]
    #         mins = i*stepd
    #     else:
    #         highs = ixd[0]
    #         mins = ixd[0] - iyd[0]
    #     ixloc = imgx_dat[mins:highs][:]
    #     print(ixloc.shape)
    #     img_product = fftpack.fft2(ixloc) * fftpack.fft2(imgy_dat).conj()
    #     t_prod = fftpack.fft(np.transpose(ixloc))*fftpack.fft2(np.transpose(imgy_dat)).conj()
    #     cc_tprod = fftpack.ifft2(t_prod)
    #     cc_image = fftpack.ifft2(img_product)
    #     cc_image.shape
    #     print(np.argmax(cc_image), np.argmax(cc_tprod))
    #     doplt(ixloc)
    #     doplt(cc_image)
    # if xarea > yarea:
    #     if ixd[0] > iyd[0]:
    #         if ixd[1] == iyd[1]:
    #             pass
    #         else:
    #             pass
    #     elif ixd[1] > iyd[1]:
    #         pass
#     img_product = fftpack.fft2(imgx_dat) * fftpack.fft2(imgy_dat).conj()
#     cc_image = fftpack.ifft2(img_product)
#     cc_image.shape()
#     image_product = np.fft.fft2(image) * np.fft.fft2(offset_image).conj()



# timgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/Testing/"
# timg1 = timgpath + 'jesusc8.jpg'
# timg2 = timgpath + 'jc8slice8.jpg'
# timg3 = timgpath + 'jc8slice8cut.jpg'
#
# arg_lt = [timg1, timg2, timg3]
# arg_lt = [timg1, timg2]
# narg = [do_imgfft2(ix) for ix in arg_lt]
# narg = [do_imzncc(ix) for ix in arg_lt]
# barg = [do_znccfft2(ix) for ix in arg_lt]
# acombs = list(itr.combinations(list(range(0,len(arg_lt))),2))
# for ac in acombs:
#     fft2_croscor(narg[ac[0]],narg[ac[1]])


# for xdat in narg:
#     xc2d = cor2d(xdat, xdat, mode='same')
#     print(xc2d.max())
#     xdot = np.dot(xdat, xdat.T)
#     print(np.average(np.abs(xdot)))
#     doplt(xdat)

# for bdat in barg:
#     xc2d = cor2d(bdat, bdat, mode='same')
#     print(xc2d.max())
#     doplt(bdat)

def get_imgdata(imgx, imgy):

    temp_list = []

    reimgx = imread(imgx).astype(float)
    reimgy = imread(imgy).astype(float)

    hx, wx, chx = np.shape(reimgx)
    print("Height: %s" % str(hx))
    print("Width: %s" % str(wx))
    print("BPP : %s " % str(chx))

    hy, wy, chy = np.shape(reimgy)
    print("Height: %s" % str(hy))
    print("Width: %s" % str(wy))
    print("BPP : %s " % str(chy))


    # gscalex = sp.inner(reimgx, [299, 587, 114]) / 1000.0
    gscalex = img_grayscale(reimgx)
    gscaley = img_grayscale(reimgy)
    ixfft = np.fft.fft2(gscalex)
    iyfft = np.fft.fft2(gscaley)

    exit()
    # ixy_cor = ixfft * (iyfft.conj())
    cc_image = np.fft.fftshift(np.fft.ifft2(ixy_cor))
    print("CCorr n x n", ixy_cor.shape)
    print("Cross Corr deg val: ", average(ixy_cor.real))
    print("Imgx FFT: ", sys.getsizeof(ixfft))
    print("Imgy FFT: ", sys.getsizeof(iyfft))


def imgproc_fft(img_data):

    pass

def img_grayscale(imgarray):
    # Do grayscale if img array is 3 dim
    if len(imgarray.shape) == 3:
        return average(imgarray, -1)
    else:
        return imgarray


def realgray_shift(imgarray):

    if len(imgarray.shape) == 3:
        return sp.inner(imgarray, [299, 587, 114]) / 1000.0
    else:
        return imgarray
