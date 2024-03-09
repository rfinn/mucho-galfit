#!/usr/bin/env python


"""

python ~/github/mucho-galfit/remove_from_mask.py maskfilename 143

this will set the pixels with values equal to 143 to zero

"""


from astropy.io import fits
import sys
import os

mfile = sys.argv[1]

value = int(sys.argv[2])


hdu = fits.open(mfile)

hdu[0].data[hdu[0].data == value] = 0

# save the prior version of the mask
os.rename(mfile,'old_'+mfile)

# write out updated mask
hdu.writeto(mfile,overwrite=True)

