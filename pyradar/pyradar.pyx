"""
Pyrex extension classes used by 'radar.py'.
"""

cdef extern from "pyradar.h":
	int radar_run_from_files(char *, char *, char *, char *, unsigned int)
	void radar_setLogLevel(int)
    
def run_from_files(filename_sequence, filename_mali,
                   filename_lfasta, filename_lfasta2,
                   random_seed):
	return radar_run_from_files(filename_sequence, filename_mali,
                                    filename_lfasta, filename_lfasta2,
                                    random_seed)
    
def setLogLevel(v):
	"""set the logging level."""
	radar_setLogLevel(v)

     
