#!/usr/bin/env python

import os, glob, sys, subprocess
import ez_setup
ez_setup.use_setuptools()

from distutils.core import setup
from distutils.extension import Extension

## check if cython exists
try:
    from Cython.Distutils import build_ext
except:
    print """Cython not installed. Please install cython from
    http://cython.org.
"""
    sys.exit (1)

## check if lfasta exitst
retcode = subprocess.call( "which lfasta",
                           shell=True)
if retcode != 0:
    print """WARNING: could not find lfasta in your path.

Please download the source code of W.R. Pearson's fasta package and build
lfasta using the following commands:

mkdir fasta2
cd fasta2
wget http://faculty.virginia.edu/wrpearson/fasta/fasta2/fasta2.shar.Z
gunzip fasta2.shar.Z
sh fasta2.shar
make lfasta

Copy the lfasta executable somewhere into your path.
"""

c_sources = [ x for x in glob.glob( "src/*.c" ) if x not in "src/main.c" ]
                   
# build the interface
pyradar = Extension(
    "pyradar",                   # name of extension
    [ "pyradar/pyradar.pyx", ] + c_sources, 
      library_dirs=[],
      libraries=[],              
      language="c",               
    )

setup(name='Radar',
      version='1.1.3',
      description='RADAR - Rapid Automatic Detection and Alignment of Repeats',
      author='Andreas Heger',
      author_email='andreas.heger@helsinki.fi',
      url='http://wwwfgu.anat.ox.ac.uk/~andreas',
      package_dir = {}, 
      packages = [], 
      scripts=['scripts/radar.py',],
      ext_modules=[ pyradar ],
      cmdclass = {'build_ext': build_ext}
     )
    