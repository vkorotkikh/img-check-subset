#!/usr/bin/python
# ******************************************************************************
# Name:    Image Subset Checker
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import sys
import scipy as sp
import numpy as np
from PIL import Image
from scipy.misc import imread
from scipy.linalg import norm
from scipy import sum, average



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

    ixy_cor = ixfft * (iyfft.conj())
    cc_image = np.fft.fftshift(np.fft.ifft2(ixy_cor))
    print("CCorr n x n", ixy_cor.shape)
    print("Cross Corr deg val: ", average(ixy_cor.real))
    print("Imgx FFT: ", sys.getsizeof(ixfft))
    print("Imgy FFT: ", sys.getsizeof(iyfft))


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
