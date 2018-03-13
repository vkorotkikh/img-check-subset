#!/usr/bin/python
# ******************************************************************************
# Name:    Image Subset Checker
# Author:  Vadim Korotkikh
# Date:    January 2018
# Info:    Program takes in two images as terminal/console arguments during exe
# and checks if either one of the images is related to another
# ******************************************************************************

import time
import logging, inspect

import numpy as np
import scipy as sp
import itertools as itr
# from pprint import pprint
from scipy import fftpack, ndimage
from scipy.misc import imread
from scipy.signal.signaltools import correlate2d as cor2d

# Set logging basicConfig level
logging.basicConfig(level=logging.INFO)
logger_ip = logging.getLogger(__name__)

# Logger handler for output
handler = logging.FileHandler('fx_img_process.log')
handler.setLevel(logging.INFO)

# Format output for logger
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)

# Add handler file to logger
logger_ip.addHandler(handler)

#>******************************************************************************
def do_grayscale(imgarr):
    """ Simplistic grayscale operation on a imgdata array """
    if len(imgarr.shape) == 3:
        return sp.average(imgarr, -1)
    else:
        return imgarr

#>******************************************************************************
def realgray_shift(imgarray):
    if len(imgarray.shape) == 3:
        return sp.inner(imgarray, [299, 587, 114]) / 1000.0
    else:
        return imgarray

#>******************************************************************************
def do_imgznorm(ipath):
    """ imread image data and perform normalization """
    imgdata = ndimage.imread(ipath, flatten=True)
    zndata = (imgdata - imgdata.mean())/ imgdata.std()
    return zndata

#>******************************************************************************
def do_imgflipnorm(ipath):
    imgdata = ndimage.imread(ipath, flatten=True)
    dataflip = np.flipud(imgdata)
    zndata = (dataflip - dataflip.mean())/dataflip.std()
    return zndata

#>******************************************************************************
def fft2_crosscorr(imgx, imgy):
    """ Function for approximating the area(s) where subimage could be located.
    Implements FFT crosscorrelation function thats calculated over each cell of
    larger image. The larger image is split into n cells that are the exact size
    of the subimage. Each cell overlaps with 1 other cells horizontally and 1
    other cells vertically.

    Assume imgx_dat is always going to be the bigger image, area wise
    imgx - filepath for larger image
    imgy - filepath for smaller image """
    '''format numpy print float output '''

    logger_ip.info('Running %s for two images' % (inspect.stack()[0][3]))
    np.set_printoptions(precision=2, linewidth=90, suppress=True)

    imgx_dat = do_imgznorm(imgx)
    imgy_dat = do_imgznorm(imgy)
    imgy_flp = do_imgflipnorm(imgy)
    ixd = imgx_dat.shape
    iyd = imgy_dat.shape
    xht = ixd[0] # imgx height/rows
    xwd = ixd[1] # imgx width/colums
    yht = iyd[0] # imgy height/rows
    ywd = iyd[1] # imgy width/columns

    imgx_area = int(xht*xwd)
    locmax_lt = []  # For storing crosscorr values
    locflip_lt = [] # For storing flip crosscorr values
    # if xht > yht and xwd == ywd: # loop over rows
    #     loc_lt = oned_slide_xcor(imgx_dat, imgy_dat, xht, yht, xwd, 'rows')
    #     # for xloc in loc_lt:
    #     #     print(xloc)
    # elif xwd > ywd and xht == yht: # loop over columns
    #     loc_lt = oned_slide_xcor(imgx_dat, imgy_dat, xwd, ywd, xht, 'cols')
    #     # for xloc in loc_lt:
    #     #     print(xloc)
    # elif xht > yht and xwd > ywd:
    ''' mulfact - Specifies cell overlap. 2 means each cell is covered by
    1/2 of next cell horizontally and 1/2 of next vertical cell.    '''
    mulfact = 2
    divnumh = mulfact*(xht//yht) # Make stepsize 1/2 subimage size in that vec
    divnumw = mulfact*(xwd//ywd)
    stepht = int(xht/divnumh)
    stepwd = int(xwd/divnumw)
    upperht, lowerht = 0, 0
    upperwd, lowerwd = 0, 0
    rngdivnumw = range(0, divnumw)

    rcindices_lt = []
    ''' precalculated fft2 (imgy data) and conjugated. Saves calc time '''
    imgy_fftc = fftpack.fft2(imgy_dat).conj()
    imgyfl_fftc = fftpack.fft2(imgy_flp).conj()

    ''' Calculate FFT crosscorrelation for each cell of main image '''
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
        temp_ihdat = imgx_dat[lowerht:upperht, :]
        for iw in rngdivnumw:
            if (iw*stepwd + ywd) < xwd:
                upperwd = iw*stepwd + ywd
                lowerwd = iw*stepwd
            else:
                upperwd = xwd
                lowerwd = xwd - ywd
            tm_rcinds.append(((lowerht, upperht), (lowerwd, upperwd)))
            # temp_dat = imgx_dat[lowerht:upperht, lowerwd:upperwd]
            temp_dat = temp_ihdat[:,lowerwd:upperwd]
            ''' try local normalization for temp image of large image '''
            # locmax = twod_slide_xcorr(temp_dat, imgy_dat)
            locmax = twod_slide_xcorrfast(temp_dat, imgy_fftc)
            locmaxfl = twod_slide_xcorrfast(temp_dat, imgyfl_fftc)
            orig_lt.append(locmax)
            flip_lt.append(locmaxfl)
        rcindices_lt.append(tm_rcinds)
        locmax_lt.append(orig_lt)
        locflip_lt.append(flip_lt)
    locmax_arr = np.asarray(locmax_lt)
    locflp_arr = np.asarray(locflip_lt)
    ''' Normalize the result cross correlation matrix for original subim '''
    locmax_normed = (locmax_arr - locmax_arr.mean()) / (np.std(locmax_arr)/2)
    ''' Normalize the result cross correlation matrix for flip subimage '''
    locmflip_norm = (locflp_arr - locflp_arr.mean()) / (np.std(locflp_arr)/2)

    # do_print_xcorrstats(locmax_normed)
    # do_print_xcorrstats(locmflip_norm)
    nmax = 4
    if imgx_area < (1680*1050):
        pass
    elif imgx_area > (1680*1050) and imgx_area < (1920*1080):
        nmax = 8
    elif imgx_area > (1920*1080):
        nmax = 13
    xcorr_sortres_lt = xcorr_result_sort(locmax_normed, rcindices_lt, nmax)
    xcorr_sortresfl_lt = xcorr_result_sort(locmflip_norm, rcindices_lt, nmax)

    loc_results={'Orig':[], 'Flip':[]}
    for xsort in xcorr_sortres_lt:
        maxvloc, rowstup, colstup = xsort[0], xsort[1], xsort[2]
        minb, ming, minr = do_focusarea_acq(imgx, imgy, maxvloc, (rowstup, colstup))
        loc_results['Orig'].append((int(minb[0]+ming[0]+minr[0]),minb, ming,minr))

    # print("Testing xcorr_sortresfl_lt")
    for flsort in xcorr_sortresfl_lt:
        maxvloc, rowstup, colstup = flsort[0], flsort[1], flsort[2]
        minb, ming, minr = do_focusarea_acq(imgx, imgy, maxvloc, (rowstup, colstup), "flip")
        loc_results['Flip'].append((int(minb[0]+ming[0]+minr[0]),minb,ming,minr))

    resprocessed_tup, rowcol_final, stat_str = do_indexloc_preprocess(loc_results)
    return resprocessed_tup, rowcol_final, stat_str



#>******************************************************************************
# def nres_slicing(lmax_array, rowcol_inds):
def xcorr_result_sort(lmax_array, rowcol_inds, nmax=5):
    """ Using the result values to cut up the original image into individual
    elements for more indepth subimage localization
    lmax_array - numpy.ndarray(nxn) - Contains the normalized maximum cross-correlation
    values computed for each subcell of main image
    rowcol_inds - List containing tuples of (height start, height end),(width start, width end)
    nmax - int value. Specifies the maximum number of 'local' maximum values to
    look for in the lmax_array. The code checks each found max value for proximity
    to already found local max values and if it's within 1 row or column away
    ie euclidean distance 1 then it is discarded    """
    rows, cols = lmax_array.shape
    # totavrg = np.average(lmax_array)
    logger_ip.info('Running %s for two images' % (inspect.stack()[0][3]))
    lmax_stdev = np.std(lmax_array)
    abmax_val = np.amax(lmax_array)
    maxval_lt = []
    for iv in range(0, nmax):
        tmp_maxval = np.partition(lmax_array.flatten(), -2)[-iv-1]
        if not maxval_lt:
            ind_tup = np.unravel_index(np.argmax(lmax_array, axis=None), (rows,cols))
            maxval_lt.append((tmp_maxval, ind_tup))
        elif maxval_lt:
            ''' set the criteria for saving a max value from the local max val
            array and and it's lmax_array cell indices location for l processing

            Minimum - tmp_maxval > lmax_stdev (std deviation of all local max values)
            Maximum - tmp_maxval < (maxval_lt[0][0] (maximum value of lmax_array)
            - lmax_stdev
            '''
            if tmp_maxval >= (maxval_lt[0][0] - (lmax_stdev)) and tmp_maxval > lmax_stdev:
                val_row, val_col = np.where(lmax_array==tmp_maxval)
                prow_rng = range(val_row[0]-1, val_row[0]+2)
                pcol_rng = range(val_col[0]-1, val_col[0]+2)
                row_lt, col_lt = [], []
                for exv in maxval_lt:
                    # erow, ecol = exv[1][0], exv[1][1]
                    row_lt.append(exv[1][0]) # row index value (0 to maxrow)
                    col_lt.append(exv[1][1]) # col index value (0 to maxcol)
                if any(ir in prow_rng for ir in row_lt) and any(ic in pcol_rng for ic in col_lt):
                    continue
                else:
                    maxval_lt.append((tmp_maxval, (val_row[0], val_col[0])))
            elif tmp_maxval >= (maxval_lt[0][0] - (3*lmax_stdev)) and tmp_maxval > lmax_stdev:
                ''' Actually this logic here is not needed, since the only difference
                between this and the if statement is 3*lmax_stdev   '''
                val_row, val_col = np.where(lmax_array==tmp_maxval)
                prow_rng = range(val_row[0]-1, val_row[0]+2)
                pcol_rng = range(val_col[0]-1, val_col[0]+2)
                row_lt, col_lt = [], []
                for exv in maxval_lt:
                    # erow, ecol = exv[1][0], exv[1][1]
                    row_lt.append(exv[1][0]) # row index value (0 to maxrow)
                    col_lt.append(exv[1][1]) # col index value (0 to maxcol)
                if any(ir in prow_rng for ir in row_lt) and any(ic in pcol_rng for ic in col_lt):
                    continue
                else:
                    maxval_lt.append((tmp_maxval, (val_row[0], val_col[0])))
    for vals in maxval_lt:
        logger_ip.info(vals)
    for vrow in lmax_array:
        logger_ip.info(vrow)
    # absmax_val = np.partition(lmax_array.flatten(), -2)[-2]
    # abthrd_val = np.partition(lmax_array.flatten(), -2)[-3]
    amax_ind = np.unravel_index(np.argmax(lmax_array, axis=None), (rows,cols))
    # print("Maximum value", abmax_val, "Type lmax_arra:", type(lmax_array))
    gslice_lt = []
    for mvtup in maxval_lt:
        gslice_res = xcorr_resgrid_slicing(lmax_array, rowcol_inds, mvtup)
        gslice_lt.append(gslice_res)
    return gslice_lt


#>******************************************************************************
def xcorr_resgrid_slicing(lmax_array, rowcol_inds, mvaltup):
    """ Using the found top 2-3 values slice up the result matrix grid for each
    one, mapping out the pixel size of the most probable subimage location for
    each value
    lmax_array - numpy.ndarray(nxn) - Contains the normalized maximum cross-correlation
    values computed for each subcell of main image
    rowcol_inds - List containing tuples of (height start, h. end), (width start, wid. end)
    mvaltup - tuple with max value and it's (row, col) indices in lmax_array
    """
    logger_ip.info('Running %s for two images' % (inspect.stack()[0][3]))
    rows, cols = lmax_array.shape
    focus_rcinds = []
    tmax_val = mvaltup[0]
    amax_ind = mvaltup[1]
    lmax_stdev = np.std(lmax_array)
    rcind = rowcol_inds[:]
    actmval_pixel = rcind[amax_ind[0]][amax_ind[1]]

    ind_minr, ind_maxr, ind_minc, ind_maxc = 0, 0, 0, 0
    ''' Getting the lmax_array matrix indices that contain good vals '''
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
    for xr in range(ind_minr, ind_maxr):
        tmp_row = []
        for xc in range(ind_minc, ind_maxc):
            # if xr==amax_ind[0] and xc==amax_ind[1]:
            #     continue
            # 5% error range from basic physics lab
            if lmax_array[xr, xc] >= (lmax_stdev-lmax_stdev*0.05):
                tmp_row.append(rcind[xr][xc])
        if tmp_row:
            focus_rcinds.append(tmp_row)
    ''' building total dimensions of focus area where subimage is most likely
        to be located '''
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
    logger_ip.info('Focus Area: %i : %i , %i : %i ' % (rowst, rowend, colst, colend))
    return actmval_pixel, (rowst, rowend), (colst, colend)

#>******************************************************************************
def do_focusarea_acq(imgx, imgy, maxinds, imgx_bounds, orient=""):
    """ Final search algorithm over predetermined area where subimage is most
    likely to be located    """
    np.set_printoptions(precision=2, linewidth=90, suppress=True)
    rowadj = 0
    coladj = 0
    # print("Imgx focus area", imgx_bounds, "Maxinds", maxinds)
    imgx_dat = ndimage.imread(imgx)
    imgy_dat = ndimage.imread(imgy)
    if orient=="flip":
        imgy_dat = np.flipud(imgy_dat)
    xrow, xcol, rgbx = imgx_dat.shape
    yrow, ycol, rgby = imgy_dat.shape

    roweq = (xrow==yrow)
    coleq = (xcol==ycol)
    rowst, rowend = imgx_bounds[0]
    colst, colend = imgx_bounds[1]

    if orient=="flip":
        imgy_dat = np.flipud(imgy_dat)

    ''' cut out the smaller focus area from the large image '''
    imgx_spec = imgx_dat[rowst:(rowend+rowadj), colst:(colend+coladj)]
    logger_ip.info('imgx_spec area - [ %i:%i , %i:%i ]' % (rowst, (rowend+rowadj), colst, (colend+coladj)))
    ''' split imagey rgb data into seperate arrays '''
    imgyb, imgyg, imgyr = imgy_dat[:, :, 0], imgy_dat[:, :, 1], imgy_dat[:, :, 2] # For RGB image

    row_range, col_range = [], []
    row_dim = imgx_spec.shape[0] - imgy_dat.shape[0]
    col_dim = imgx_spec.shape[1] - imgy_dat.shape[1]
    if row_dim > 0:
        row_range = list(range(0,(imgx_spec.shape[0] - imgy_dat.shape[0]),2))
    elif row_dim==0:
        # row_range = range(0, 2, 1)
        row_range = range(0,1)
    if col_dim > 0:
        col_range = list(range(0,(imgx_spec.shape[1] - imgy_dat.shape[1]),2))
    else:
        # col_range = range(0, 2, 1)
        col_range = range(0,1)
    mse_ablue = []
    mse_ared = []
    mse_agre = []
    tmpr_end, tmpc_end = 0, 0
    for xr in row_range:
        tmpr_end = yrow + xr
        mse_tb = []
        mse_tr = []
        mse_tg = []
        hist_tr = []
        for xc in col_range:
            tmpc_end = ycol + xc
            loc_ispecx = imgx_spec[xr:tmpr_end, xc:tmpc_end]
            # if (yrow, ycol) == loc_ispecx.shape[0:2]:
            #     pass
            # else:
            #     continue
            mse_blue = do_mse(loc_ispecx[:,:,0], imgyb)
            mse_gre = do_mse(loc_ispecx[:,:,1], imgyg)
            mse_red = do_mse(loc_ispecx[:,:,2], imgyr)
            mse_tb.append(mse_blue)
            mse_tr.append(mse_red)
            mse_tg.append(mse_gre)
        mse_ablue.append(mse_tb)
        mse_agre.append(mse_tg)
        mse_ared.append(mse_tr)
    npar_blue = np.asarray(mse_ablue)
    npar_green = np.asarray(mse_agre)
    npar_red = np.asarray(mse_ared)
    # mxrows = maxinds[0]
    # mxcols = maxinds[1]
    # adjxrows = (mxrows[0] - rowst, mxrows[1] - rowst)
    # adjxcols = (mxcols[0] - colst, mxcols[1] - colst)

    mse_allshp = npar_blue.shape
    nblueixs = np.unravel_index(np.argmin(npar_blue, axis=None), mse_allshp)
    radjbr, radjbc = rowst+ 2*nblueixs[0],  colst + 2*nblueixs[1]
    ngreenixs = np.unravel_index(np.argmin(npar_green, axis=None), npar_green.shape)
    radjgr, radjgc = rowst + 2*ngreenixs[0], colst + 2*ngreenixs[1]
    nredixs = np.unravel_index(np.argmin(npar_red, axis=None), npar_red.shape)
    radjr, radjc = rowst + 2*nredixs[0], colst + 2*nredixs[1]
    nminblue, nmingreen, nminred = np.min(npar_blue), np.min(npar_green), np.min(npar_red)
    # print(nminblue, nblueixs, (radjbr, radjbc))
    # print(nmingreen, ngreenixs, (radjgr, radjgc))
    # print(nminred, nredixs, (radjr, radjc))

    return (nminblue, nblueixs, (radjbr, radjbc)), (nmingreen, ngreenixs, (radjgr, radjgc)), (nminred, nredixs, (radjr, radjc))

#>******************************************************************************
def do_indexloc_preprocess(locresults_dict):
    """ Run results processing over all possible image orinetations
    Currently just original and flipped     """
    reskeys = locresults_dict.keys()

    subimage_stat = ""
    row_ind, col_ind = 0, 0
    comb_minvals = []
    minv_dict = {}
    for rkey in reskeys:
        rlist = locresults_dict[rkey]
        sumvals = [x[0] for x in rlist]
        indxmin = np.argmin(sumvals)
        minv_dict[rkey] = rlist[indxmin]
        comb_minvals.append((min(sumvals), rkey))
    finmin = min(minv_dict.items(), key=lambda x: x[1][0])
    pixel_inds = [x[2] for x in finmin[1][1:]]

    if all(x[0] == 1 for x in pixel_inds) and (x[1] == 1 for x in pixel_inds):
        row_ind = int(np.median([xi[0] for xi in pixel_inds]))
        col_ind = int(np.median([xi[1] for xi in pixel_inds]))
        subimage_stat = "True"
    elif pixel_inds[0][0] == pixel_inds[1][0] and pixel_inds[0][1] == pixel_inds[1][1]:
        row_ind = int(np.median([xi[0] for xi in pixel_inds]))
        col_ind = int(np.median([xi[1] for xi in pixel_inds]))
        subimage_stat = "True"
    elif pixel_inds[1][0] == pixel_inds[2][0] and pixel_inds[1][1] == pixel_inds[2][1]:
        row_ind = int(np.median([xi[0] for xi in pixel_inds]))
        col_ind = int(np.median([xi[1] for xi in pixel_inds]))
        subimage_stat = "True"
    elif pixel_inds[0][0] == pixel_inds[2][0] and pixel_inds[0][1] == pixel_inds[2][1]:
        row_ind = int(np.median([xi[0] for xi in pixel_inds]))
        col_ind = int(np.median([xi[1] for xi in pixel_inds]))
        subimage_stat = "True"
    elif pixel_inds[0] != pixel_inds[1] != pixel_inds[2]:
        # print("Not a subimage")
        subimage_stat = "False"
    else:
        pass
    # print(pixel_inds)

    return finmin, (row_ind, col_ind), subimage_stat

#>******************************************************************************
def do_mse(imageA, imageB):
    # the 'Mean Squared Error' between the two images is the
    # sum of the squared difference between the two images;
    err = np.sum((imageA.astype("float") - imageB.astype("float")) ** 2)
    err /= float(imageA.shape[0] * imageA.shape[1])
    return err

#>******************************************************************************
def twod_slide_xcorr(islice, subimg):
    img_product = fftpack.fft2(islice) * fftpack.fft2(subimg).conj()
    inv_prod = fftpack.ifft2(img_product)
    return np.amax(inv_prod.real)

#>******************************************************************************
def twod_slide_xcorrfast(islice, simgfcon):
    """ Execute cross correlation """
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
def do_print_xcorrstats(locmax_normvals):
    # print(np.mean(locmax_normvals), np.average(locmax_normvals))
    print("")
    print("Total Avg: ", np.average(locmax_normvals), "Standard Deviation", np.std(locmax_normvals))
    for xarr in locmax_normvals:
        print(xarr, "", np.average(xarr))

#>******************************************************************************
def get_npstats(np_onedarray):
    """ Simple function to get and return avg, maxvalue and standard deviation
    values for a given np array """
    arr_avg = np.average(np_onedarray)
    arr_amax = np.average(np_onedarray)
    arr_stdev = np.std(np_onedarray)
    return arr_avg, arr_amax, arr_stdev
