/*
RADAR - a program to rapidly detect and align repeats in protein sequences
Copyright (C) 2000 Andreas Heger

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

this program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

##########################################################################
$Id: radar.c,v 1.2 2005/02/17 23:13:32 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository
			Fixed bug in AlignBestWindow => origdots is copied
30.9.2003       heger   General overhaul. Removed some parameters. Deduce
                        length of sequence/mali from file.
 */

#include "radar.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>

extern VERBOSITY verbose;

/* global parameters */

static char filename_ma[100];
static char filename_lfasta[100];
static char filename_lfasta2[100];
static char filename_sequence[100];
static int random_seed=0;

/* Version 7: 
   -> collect further instances of repeats by simple alignment without wrapping
   -> reduce diagonal by mapping dots for bestwindow, but not for alignment
   -> calculate register, align repeats, and then calculate register again
   -> reduce profile
   -> if you do not find a bestwindow, look at maxdiagonals best diagonals
 */

/* Version13:
   use FastaFiles2DotsSaveSpace instead of FastaFiles2DotsCorrect,
   somehow not equivalent results to 12, check later
 */

/* Version 15a: release for Rodrigo, changed output */

/*--------------start of parameter parsing -----------------------------------*/
const char * my_progname = "radar";

static void print_version() {
	printf("%s Version $Id: radar.c,v 1.2 2005/02/17 23:13:32 aheger Exp $\n", my_progname);
}    

static void usage()
{
	print_version();
	printf("Usage: %s [PARAMETERS] [OPTIONS]\n", my_progname);

	printf("-V		print version and exit.\n");
	printf("-h		print help and exit.\n");
	printf("-v #		loglevel (optional).\n");
	printf("-P #		filename of multiple alignment [%s](required).\n", filename_ma);
	printf("-Q #		filename of lfasta [%s] (required).\n", filename_lfasta);
	printf("-R #		filename of sequence [%s] (required).\n", filename_sequence);
	printf("-S #		filename of lfasta [%s] (optional).\n", filename_lfasta2);
}

void ParseArguments (int argc, char *argv[]) {

	int c;  
	extern char * optarg;

	while ((c=getopt(argc, argv, "?Vv:A:B:P:Q:R:S:")) != EOF) {
		switch(c) {
		case 'V':
			print_version(); exit(EXIT_SUCCESS);
		case 'h':
		case '?':
			usage(); exit(EXIT_SUCCESS);
		case 'v':
			verbose = atoi(optarg); break;
		case 'P':	    /* filename sequence1 */
			strcpy(filename_ma, optarg); break;
		case 'Q':    	    /* filename sequence2 */
			strcpy(filename_lfasta, optarg); break;
		case 'R':    	    /* filename sequence3 */
			strcpy(filename_sequence, optarg); break;
		case 'S':    	    /* filename sequence3 */
			strcpy(filename_lfasta2, optarg); break;
		}
	}

	// set pointers to end of options
	(argc)-=optind;
	(argv)+=optind;

	if (argc > 0 || 
			strlen(filename_ma) == 0 ||
			strlen(filename_lfasta) == 0 ||
			strlen(filename_sequence) == 0) {
		usage();
		exit(EXIT_FAILURE);
	}
}


/*-------------------------------------------------------------------*/	    
/* subroutine does the iteration */
/*--------------------------------------------------------------------------*/
int main(argc, argv)
int argc; char *argv[];
{
	ParseArguments(argc, argv);
	return radar_run_from_files(filename_sequence,
				    filename_ma,
				    filename_lfasta,
				    filename_lfasta2,
				    random_seed);
}
