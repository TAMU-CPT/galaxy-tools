#!/usr/bin/env python
import h5py
import numpy
import sys

OUT_FILE_PATH = sys.argv[1]
PARENT_BASH5 = sys.argv[2]
CHILDREN = sys.argv[3:]


import os
import shutil
if not os.path.exists(OUT_FILE_PATH):
    os.makedirs(OUT_FILE_PATH)

NEW_BAS_H5 = os.path.join(OUT_FILE_PATH, 'parent.bas.h5')
shutil.copy(PARENT_BASH5, NEW_BAS_H5)
for child in CHILDREN:
    shutil.copy(child, os.path.join(OUT_FILE_PATH, child.rsplit('/', 1)[1]) +
                '.bax.h5')

f = h5py.File(NEW_BAS_H5, 'r+')
# Remove the original
try:
    del f['MultiPart/Parts']
except:
    pass

dt = h5py.special_dtype(vlen=bytes)
# Re-add the dataset with the names of whichever children are passed in
dset = f.create_dataset('MultiPart/Parts',
                        data=numpy.array([x.rsplit('/', 1)[1] + '.bax.h5' for x
                                          in CHILDREN]), dtype=dt)
# Add the description though they probably don't need it/won't check it
dset.attrs.create('Description', numpy.array(['File part names']), dtype=dt)
