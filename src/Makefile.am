# Copyright (C) 2000 Andreas Heger <heger@ebi.ac.uk>
#  
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without 
# modifications, as long as this notice is preserved.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


INCLUDES = -I$(top_srcdir)/src 

EXTRA_DIST=radar.h align.h align_io.h align_lib.h toolbox.h align_tools.h

bin_PROGRAMS = radar

radar_SOURCES = radar.c register.c repeatstatistics.c \
	dotalign.c dotalign_tube.c fasta2dot.c \
	toolbox.c align_io.c align_tools.c counts2profile.c \
	getmaxdiag.c lps_align.c gps_align.c 
