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
def do_imgznorm(ipath):
    """ imread image data and perform normalization """
    imgdata = ndimage.imread(ipath, flatten=True)
    zndata = (imgdata - imgdata.mean())/ imgdata.std()
    # zndata = (imgdata - imgdata.mean())
    return zndata

#>******************************************************************************
def do_imgflipnorm(ipath):
    imgdata = ndimage.imread(ipath, flatten=True)
    dataflip = np.flipud(imgdata)
    zndata = (dataflip - dataflip.mean())/dataflip.std()
    return zndata

#>******************************************************************************
def fft2_crosscorr(imgx, imgy):
    """ Assume imgx_dat is always going to be the bigger image, area wise
    imgx - filepath for larger imagex
    imgy - filepath for smaller imagey """
    '''format float output from calc '''
    np.set_printoptions(precision=2, linewidth=90, suppress=True)

    imgx_dat = do_imgznorm(imgx)
    imgy_dat = do_imgznorm(imgy)
    imgy_flp = do_imgflipnorm(imgy)
    imgx_dim = imgx_dat.shape
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
    locflip_lt = [] # For storring flip comparison values
    if xht > yht and xwd == ywd: # loop over rows
        loc_lt = oned_slide_xcor(imgx_dat, imgy_dat, xht, yht, xwd, 'rows')
        for xloc in loc_lt:
            print(xloc)
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

        imgy_fftc = fftpack.fft2(imgy_dat).conj()
        imgyfl_fftc = fftpack.fft2(imgy_flp).conj()
        print("Ht",xht,yht," Wd", xwd, ywd)
        print("Stepht: ", stepht, "Stephwd: ", stepwd)
        """ Build the index - pixel loc matrix """
        rcindices_lt = []
        for ih in range(0, divnumh):
            tm_rcinds = []
            orig_lt = []
            flip_lt = []
            if (ih*stepht + yht) < xht:
                upperht = ih*stepht + yht
                lowerht = ih*stepht
            else:
                upperht = xht
                lowerht = xht - yht
            # rcindices_lt.append((lowerht, upperht))
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
                ''' try local normalization for temp image of large image '''
                # locmax = twod_slide_xcorr(temp_dat, imgy_dat)
                locmax = twod_slide_xcorrfast(temp_dat, imgy_fftc)
                # locmaxfl = twod_slide_xcorr(temp_dat, imgy_flp)
                locmaxfl = twod_slide_xcorrfast(temp_dat, imgyfl_fftc)
                # print("Indices", rowstr, colstr, " LocMax: ", locmax)
                orig_lt.append(locmax)
                flip_lt.append(locmaxfl)
            rcindices_lt.append(tm_rcinds)
            locmax_lt.append(orig_lt)
            locflip_lt.append(flip_lt)

        locmax_arr = np.asarray(locmax_lt)
        locflp_arr = np.asarray(locflip_lt)
        print(np.mean(locmax_arr), np.average(locmax_arr))
        ''' Normalize the result cross correlation matrix for original subim '''
        locmax_normed = (locmax_arr - locmax_arr.mean()) / (np.std(locmax_arr)/2)
        ''' Normalize the result cross correlation matrix for flip subimage '''
        locmflip_norm = (locflp_arr - locflp_arr.mean()) / (np.std(locflp_arr)/2)
        print(np.mean(locmax_normed), np.average(locmax_normed))
        print("")
        print("Total Avg: ", np.average(locmax_normed), "Stdev", np.std(locmax_normed))
        for xarr in locmax_normed:
            print(xarr, "", np.average(xarr))
        print("")
        print("Total Avg: ", np.average(locmflip_norm), "Stdev", np.std(locmflip_norm))
        for xfl in locmflip_norm:
            print(xfl,"", np.average(xfl))
        maxvloc, rowstup, colstup = nres_slicing(locmax_normed, rcindices_lt)
        maxfloc, rowsflp, colsflip = nres_slicing(locmflip_norm, rcindices_lt)

        minrgb = do_focusarea_acq(imgx, imgy, maxvloc, (rowstup, colstup))
        minrgbflip = do_focusarea_acq(imgx, imgy, maxfloc, (rowsflp, colsflip), "flip")
        # return rowstup, colstup



#>******************************************************************************
def nres_slicing(lmax_array, rowcol_inds):
    """ Using the result values to cut up the original image into individual
    elements for more indepth subimage localization
    lmax_array - Contains the normalized maximum coefficient values computed for
    each subcell of main image
    rowcol_inds - List containing tuples of (height start, height end), (width"""
    rcind = rowcol_inds[:]
    rows, cols = lmax_array.shape
    totavrg = np.average(lmax_array)
    totstdev = np.std(lmax_array)
    amax_arr = np.amax(lmax_array)
    # am_ind = np.argmax(lmax_array)
    nlarge = np.partition(lmax_array.flatten(), -2)[-2]
    amax_ind = np.unravel_index(np.argmax(lmax_array, axis=None), (rows,cols))
    actmax_pixel = rcind[amax_ind[0]][amax_ind[1]]
    print("2nd", nlarge)
    print("Amax index:", amax_ind)
    print(rcind[amax_ind[0]][amax_ind[1]])
    # amax_stdev = amax_arr//totstdev
    focus_rcinds = []
    # nw_inds = [0, 0]
    # ew_inds = [0, 0]
    # sw_inds = [0, 0]
    # se_inds = [0, 0]
    ind_minr, ind_maxr, ind_minc, ind_maxc = 0, 0, 0, 0
    if amax_ind[0] >= 1:
        ind_minr = amax_ind[0]-1
        if (amax_ind[0]+2) <= rows:
            ind_maxr = amax_ind[0]+2
        else:
            ind_maxr = amax_ind[0]+1
        # ind_minr, ind_maxr = amax_ind[0]-1, amax_ind[0]+2
    else:
        ind_minr = amax_ind[0]
        if (amax_ind[0]+2) <= rows:
            ind_maxr = amax_ind[0]+2
        else:
            ind_maxr = amax_ind[0]+1

    if amax_ind[1] >= 1:
        ind_minc = amax_ind[1]-1
        if (amax_ind[1]+2)<=cols:
            ind_maxc =  amax_ind[1]+2
        else:
            ind_maxc =  amax_ind[1]+1
    else:
        ind_minc = amax_ind[1]
        if (amax_ind[1]+2)<=cols:
            ind_maxc =  amax_ind[1]+2
        else:
            ind_maxc =  amax_ind[1]+1
    # ind_minc, ind_maxc = amax_ind[1]-1, amax_ind[1]+2
    print(ind_minr, ind_maxr)
    print(ind_minc, ind_maxc)
    for xr in range(ind_minr, ind_maxr):
        tmp_row = []
        for xc in range(ind_minc, ind_maxc):
            # if xr==amax_ind[0] and xc==amax_ind[1]:
            #     continue
            # 5% error range from basic physics lab
            if lmax_array[xr, xc] >= (totstdev-totstdev*0.05):
                tmp_row.append(rcind[xr][xc])
        if tmp_row:
            focus_rcinds.append(tmp_row)
    ''' building total dimensions of focus area where subimage is most likely
        to be located '''
    for xrc in focus_rcinds:
        print(xrc)
    rowst, rowend = 0, 0
    colst, colend = 0, 0
    # rowst = focus_rcinds[0][0][0]
    for rcx in focus_rcinds:
        # print(rcx[0], "", rcx[-1])
        if rcx[0][0][0] > 0 and rowst == 0:
            rowst = rcx[0][0][0]
        if rcx[0][1][0] > 0 and colst == 0:
            colst = rcx[0][1][0]
        if rcx[0][0][1] > rowend:
            rowend = rcx[0][0][1]
        if rcx[-1][1][1] > colend:
            colend = rcx[-1][1][1]
        # rowst, colst = rcx[0][0][0], rxc[0][1][0]
        # print(rcx)
    print(rowst, rowend)
    print(colst, colend)
    return actmax_pixel, (rowst, rowend), (colst, colend)

#>******************************************************************************
def do_focusarea_acq(imgx, imgy, maxinds, imgx_bounds, orient=""):
    """ Final search algorithm over predetermined area where subimage is most
    likely to be located    """
    np.set_printoptions(precision=2, linewidth=90, suppress=True)
    rowadj = 0
    coladj = 0
    print("Imgx focus area", imgx_bounds, "Maxinds", maxinds)
    imgx_dat = ndimage.imread(imgx)
    imgy_dat = ndimage.imread(imgy)
    if orient=="flip":
        imgy_dat = np.flipud(imgy_dat)
    yrow, ycol, rgbval = imgy_dat.shape
    rowst, rowend = imgx_bounds[0]
    colst, colend = imgx_bounds[1]
    imgx_spec = imgx_dat[rowst:(rowend+rowadj), colst:(colend+coladj)]
    mxrows = maxinds[0]
    mxcols = maxinds[1]
    adjxrows = (mxrows[0] - rowst, mxrows[1] - rowst)
    adjxcols = (mxcols[0] - colst, mxcols[1] - colst)
    print("Adj", adjxrows, adjxcols)
    print(imgx_spec.shape)
    # print(imgx_spec[0,0])
    # print(imgx_spec[0,0:5])
    imgyb, imgyg, imgyr = imgy_dat[:, :, 0], imgy_dat[:, :, 1], imgy_dat[:, :, 2] # For RGB image
    for c in range(0,3):
        tmp_im = np.zeros(imgx_spec.shape, dtype="uint8")
        tmp_imy = np.zeros(imgy_dat.shape, dtype="uint8")
        tmp_im[:,:,c] = imgx_spec[:,:,c]
        # imgy_dat[:,:,c] = (imgy_dat[:,:,c] - imgy_dat[:,:,c].mean())/np.std(imgy_dat[:,:,c])
        # imgx_spec[:,:,c] = (imgx_spec[:,:,c] - imgx_spec[:,:,c].mean())/np.std(imgx_spec[:,:,c])

    row_range, col_range = [], []
    rowrlol = imgx_spec.shape[0] - imgy_dat.shape[0]
    colrlol = imgx_spec.shape[1] - imgy_dat.shape[1]
    if rowrlol > 0:
        row_range = list(range(0,(imgx_spec.shape[0] - imgy_dat.shape[0]),2))
    else:
        row_range = range(0, 1)
    if colrlol > 0:
        col_range = list(range(0,(imgx_spec.shape[1] - imgy_dat.shape[1]),2))
    else:
        col_range = range(0, 1)
    mse_ablue = []
    mse_ared = []
    mse_agre = []
    print(rowrlol, colrlol)
    for xr in row_range:
        rowend = yrow + xr
        # if xr > 10:
        #     exit()
        mse_tb = []
        mse_tr = []
        mse_tg = []
        for xc in col_range:
            colend = ycol + xc
            loc_ispecx = imgx_spec[xr:rowend, xc:colend]
            mse_blue = do_mse(loc_ispecx[:,:,0], imgyb)
            mse_gre = do_mse(loc_ispecx[:,:,1], imgyg)
            mse_red = do_mse(loc_ispecx[:,:,2], imgyr)
            mse_tb.append(mse_blue)
            mse_tr.append(mse_red)
            mse_tg.append(mse_gre)
            # print("Imgx Var", np.var(loc_ispecx[:,:,0]), "Std", np.std(loc_ispecx[:,:,0]), "Avg", np.mean(loc_ispecx[:,:,0])) # b
            # print("Imgx Var", np.var(loc_ispecx[:,:,1]), "Std", np.std(loc_ispecx[:,:,1]), "Avg", np.mean(loc_ispecx[:,:,1])) # g
            # print("Imgx Var", np.var(loc_ispecx[:,:,2]), "Std", np.std(loc_ispecx[:,:,2]), "Avg", np.mean(loc_ispecx[:,:,2])) # r
        mse_ablue.append(mse_tb)
        mse_ared.append(mse_tr)
        mse_agre.append(mse_tg)
    npar_blue = np.asarray(mse_ablue)
    npar_green = np.asarray(mse_agre)
    npar_red = np.asarray(mse_ared)
    # for xar in npar_blue:
    #     print(xar)
    # print("")
    # for xgre in npar_green:
    #     print(xgre)
    # print("")
    # for xar in npar_red:
    #     print(xar)
    mse_allshp = npar_blue.shape
    nminblue, nmingreen, nminred = np.min(npar_blue), np.min(npar_green), np.min(npar_red)
    print(nminblue, np.unravel_index(np.argmin(npar_blue, axis=None), mse_allshp))
    print(nmingreen, np.unravel_index(np.argmin(npar_green, axis=None), mse_allshp))
    print(nminred, np.unravel_index(np.argmin(npar_red, axis=None), mse_allshp))
    return nminblue, nmingreen, nminred

#>******************************************************************************
def do_mse(imageA, imageB):
	# the 'Mean Squared Error' between the two images is the
	# sum of the squared difference between the two images;
	err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
	err /= float(imageA.shape[0] * imageA.shape[1])

	return err

#>******************************************************************************
def get_npstats(np_onedarray):
    """ Simple function to get and return avg, maxvalue and standard deviation
    values for a given np array """
    arr_avg = np.average(np_onedarray)
    arr_amax = np.average(np_onedarray)
    arr_stdev = np.std(np_onedarray)
    return arr_avg, arr_amax, arr_stdev


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
def realgray_shift(imgarray):

    if len(imgarray.shape) == 3:
        return sp.inner(imgarray, [299, 587, 114]) / 1000.0
    else:
        return imgarray
