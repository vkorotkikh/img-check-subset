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
elif pyver_tup >= (2, 7) and pyver_tup <= (3,0):
    pass
else:
    raise Exception("Minimum Python ver 3.5 or 2.7 are required")


#>******************************************************************************
def main(imgpath_x, imgpath_y):

    if checkfile(imgpath_x):
        if checkfile(imgpath_y):
            pass
        else:
            sys.exit("Image file DNE")
    else:
        sys.exit("Image file DNE")

    print("Python %s.%s" % (pyver_major, pyver_minor))

#>******************************************************************************
def checkfile(ifilepath):
    if os.path.isfile(ifilepath):
        return 1
    else:
        return 0




if __name__ == "__main__":
    try:
        main(sys.argv[1], sys.argv[2])
    except IndexError:
        sys.exit("Two Terminal arguments are required")
