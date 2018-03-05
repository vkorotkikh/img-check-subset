#!/usr/bin/python
# ******************************************************************************
# Name: Class ImgData
# Desig: Contains various functions for storing basic class info
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import os, sys
import numpy as np

class ImgData():


    def __init__(self, imagepath, imgdata):`
        self.imagepath = imagepath
        self.imgdata = imgdata
        self.area = defarg
        # self.cols =
        # self.rows =



    #>******************************************************************************
    def do_grayscale(self, imgarr):
        """ Simplistic grayscale operation on a imgdata array """
        if len(imgarr.shape) == 3:
            return sp.average(imgarr, -1)
        else:
            return imgarr

    #>******************************************************************************
    def do_imgznorm(self, ipath):
        """ imread image data and perform normalization """
        imgdata = ndimage.imread(ipath, flatten=True)
        zndata = (imgdata - imgdata.mean())/ imgdata.std()
        # zndata = (imgdata - imgdata.mean())
        return zndata

    #>******************************************************************************
    def do_imgflipnorm(self, ipath):
        imgdata = ndimage.imread(ipath, flatten=True)
        dataflip = np.flipud(imgdata)
        zndata = (dataflip - dataflip.mean())/dataflip.std()
        return zndata
