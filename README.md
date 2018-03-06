# Image - Subimage Localization Utility
This tool is used to determine whether two images are related and whether one
is a subimage of the other.

Uses a combination of FFT based normalized cross-correlation method and MSE to
perform subimage localization.

## Important details
Developed and tested using Python 3.6 using Anaconda Conda installation.
Requires numpy and scipy libraries
```
$ Requires numpy and scipy libraries
$ Tested working with Python 3.5 in MacOSX environment using Sierra 10.12.6
```

## Getting Started
To execute this tool, run the following command in Terminal
```
$ python -u rung_imgchk.py .../image1 .../image2
```
or if you have both Python 2 and 3, specify Python 2 version.
```
$ python3 -u rung_imgchk.py .../image1 .../image2
```

## Author
Vadim Korotkikh
# img-check-subset repository
