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
from scipy import fftpack, ndimage, average
from scipy.misc import imread
from scipy.signal.signaltools import correlate2d as cor2d


def do_grayscale(imgarr):
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
def fft2_croscor(imgx_dat, imgy_dat):
    """ Assume imgx_dat is always going to be the bigger image, area wise"""
    imgx_dim = imgx_dat.shape
    imgy_dim = imgy_dat.shape
    ixd, iyd = imgx_dim, imgy_dim
    xht = ixd[0] # heigth
    xwd = ixd[1] # width
    yht = iyd[0] # h
    ywd = iyd[1] # w

    divnum = 0
    if xht > yht and xwd == ywd:
        divnum = xht//yht
        stepd = 4*int(xht/divnum)
        for i in range(0, divnum):
            upperb = 0
            lowerb = 0
            if i*stepd + yht < xht:
                upperb = i*stepd + yht
                lowerb = i*stepd
            else:
                upperb = xht
                lowerb = xht - yht
    elif xwd > ywd and xht == yht:
        divum = xwd//ywd
        stepd = 4*int(xwd/divnum)
        for i in range(0, divnum):
            upperb = 0
            lowerb = 0
            if i*stepd + ywd < xwd:
                upperb = i*stepd + ywd
                lowerb = i*stepd
            else:
                upperb = xwd
                lowerb = xwd - ywd

    elif xht > yht and xwd > ywd:
        ''' this one is messy '''
        pass

    divn = 4*(ixd[0]//iyd[0])
    stepd = int(ixd[0]/divn)

    for i in range(0, divn):
        highs = 0
        if i*stepd + iyd[0] < ixd[0]:
            highs = i*stepd + iyd[0]
            mins = i*stepd
        else:
            highs = ixd[0]
            mins = ixd[0] - iyd[0]
        ixloc = imgx_dat[mins:highs][:]
        print(ixloc.shape)
        img_product = fftpack.fft2(ixloc) * fftpack.fft2(imgy_dat).conj()
        t_prod = fftpack.fft(np.transpose(ixloc))*fftpack.fft2(np.transpose(imgy_dat)).conj()
        cc_tprod = fftpack.ifft2(t_prod)
        cc_image = fftpack.ifft2(img_product)
        cc_image.shape
        print(np.argmax(cc_image), np.argmax(cc_tprod))
        doplt(ixloc)
        doplt(cc_image)
    if xarea > yarea:
        if ixd[0] > iyd[0]:
            if ixd[1] == iyd[1]:
                pass
            else:
                pass
        elif ixd[1] > iyd[1]:
            pass
#     img_product = fftpack.fft2(imgx_dat) * fftpack.fft2(imgy_dat).conj()
#     cc_image = fftpack.ifft2(img_product)
#     cc_image.shape()
#     image_product = np.fft.fft2(image) * np.fft.fft2(offset_image).conj()



timgpath = "/Users/vkorotki/Movies/Utils/img-check-subset/Testing/"
timg1 = timgpath + 'jesusc8.jpg'
timg2 = timgpath + 'jc8slice8.jpg'
timg3 = timgpath + 'jc8slice8cut.jpg'

# arg_lt = [timg1, timg2, timg3]
arg_lt = [timg1, timg2]
# narg = [do_imgfft2(ix) for ix in arg_lt]
narg = [do_imzncc(ix) for ix in arg_lt]
barg = [do_znccfft2(ix) for ix in arg_lt]
acombs = list(itr.combinations(list(range(0,len(arg_lt))),2))
for ac in acombs:
    fft2_croscor(narg[ac[0]],narg[ac[1]])


# for xdat in narg:
#     xc2d = cor2d(xdat, xdat, mode='same')
#     print(xc2d.max())
#     xdot = np.dot(xdat, xdat.T)
#     print(np.average(np.abs(xdot)))
#     doplt(xdat)
print("\n")
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
