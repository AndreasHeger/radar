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

/* gets: ma and fasta-file and sequence */

/* global exported variables - used in main.c */
VERBOSITY verbose = LL0;

char globalresult[MAXRESULTSIZE];

/* global parameters */
static int mindiag		= 15;   /* minimum distance for dots from diagonal */
static int maxlevels           	= 1;    /* max level for sampling diagonals */
static int maxdiagonals        	= 3;    /* max level for sampling diagonals */
static int diagwidth		= 1;    /* width of window for sampling of diagonals */
static int scalefactor		= 20;   /* scalefactor for determining the width of tube in tubalignment */
static int minsamplediag	= 3;    /* minimum half-width of tube for tube alignment */
static int minlrepeat           = 10;   /* minimum repeat length */
static int minndots             = 20;   /* minimum number of dots needed for alignment */
static int splitdistance        = 50;   /* residue-separation, above which two repeats are regarded as distinct */
static int shuffle_iterations   = 100;  /* number of iterations for shuffling sequences */
static CALCTYPE selectfactor    = 0.25;  /* cutoff, that decides, which repeats are outliers and which are not; 25% of mean score */
static CALCTYPE mintubescore    = 10;   /* minimum score for tube alignment */
static CALCTYPE minswscore      = 10;   /* minimum score for wrapped alignment */
static CALCTYPE cutoffscore     = 20;	/* cutoff in PruneMA and SelectRepeats for keeping repeats*/
static CALCTYPE mincollectscore = 20;   /* minimum score for collecting repeats */
static CALCTYPE mindiagscore    = 1;    /* minimum score for diagonal to be kept (10 is too high)*/   
static CALCTYPE minzscore       = 6.0; /* minimum zscore for keeping repeats , default: 6.0*/   
static CALCTYPE reducefactor    = 0.75; /* maximum length of reduced repeat, i.e. maximum 3/4 of original repeat-length */
/* static CALCTYPE increasefactor  = 1.5;  minimum number of repeats more for reduction, i.e. keep alignment, if at least 1 1/2 more than original repeats */

/* gap penalties */

/* for reducing diagonal */
static CALCTYPE rd_gapopen_i      = -4;     
static CALCTYPE rd_gapelon_i      = -0.4;   
static CALCTYPE rd_gapopen_j      = -4;     
static CALCTYPE rd_gapelon_j      = -0.4;   

/* for register, look there as well */
static CALCTYPE rg_gapopen = -4;
static CALCTYPE rg_gapelon = -0.4; 

/* for first dotalignment */
static CALCTYPE gapopen_i      = -4;       /* gap penalties for bestwindow/repeat */
static CALCTYPE gapelon_i      = -0.4; 
static CALCTYPE gapopen_j      = -4;       /* gap penalties for sequence */   
static CALCTYPE gapelon_j      = -0.4; 

/* for collecting repeats */
static CALCTYPE hp_gapopen_i      = -4;       /* gap penalties for bestwindow/repeat */
static CALCTYPE hp_gapelon_i      = -0.4; 
static CALCTYPE hp_gapopen_j      = -4;       /* gap penalties for sequence */   
static CALCTYPE hp_gapelon_j      = -0.4; 

/* for bestwindow */
static CALCTYPE bw_gapopen_i      = -4;     /*  gap penalties for profile */
static CALCTYPE bw_gapelon_i      = -0.4;   
static CALCTYPE bw_gapopen_j      = -4;    /*  gap penalties for sequence */
static CALCTYPE bw_gapelon_j      = -0.4;     

static int iteration;

#define GAPCHAR '.';


/*--------------end of parameter parsing -----------------------------------*/

void radar_setLogLevel(int loglevel)
{
	verbose = loglevel;
}

/*--------------------------------------------------------------------------*/
void DumpAlignment( ali, title )
ALI * ali;
const char * title;
{
	int i;

	printf("---------- %s - START ---------------\n", title); 
	for (i = 0; i <= ali->lastindex; i++) {
		printf("|%5i%5i%5.2f", ali->align_i[i], ali->align_j[i], ali->align_s[i]);
	}
	printf("\n---------- %s - END ----------------\n", title); 

}

void DumpArrayAlignment( ali, title )
ARRAYALI * ali;
const char *title;
{
	int i;
	printf("---------- %s - START ---------------\n", title); 
	for (i = 1; i <= ali->length; i++) {
		printf("|%5i%5i%5.2f", i, ali->ali[i], ali->scores[i]);
	}
	printf("\n---------- %s - END ----------------\n", title); 

}

/* memory has to be available for the new dots */
void AddDiagonal ( dots, profile, frequencies, lprofile, sequence, lsequence )
DOTS *dots;
PROFILECOLUMN   *profile;
FREQUENCYCOLUMN *frequencies;
int lprofile;
SEQTYPE *sequence;
int lsequence;
{
	int i, x;
	int y = dots->ndots;
	for (i = 0; i < lsequence; i++) {
		x = i+1;
		dots->row[y]   = x;
		dots->col[y]   = x;
		dots->score[y] = ProfileScore( profile[x], frequencies[x], profile[x], frequencies[x]);
		y++;
	}
	dots->ndotsdiagonal = y;
}

SEQTYPE *AllocateSequenceMemory( size ) 
int size;
{
	SEQTYPE * seq = Malloc( sizeof( SEQTYPE) * (size + 1));
	return seq;
}

void FreeSequenceMemory( seq ) 
SEQTYPE * seq;
{
	free( seq );
}


ALI *AllocateAliMemory ( size)
int size;
{
	ALI * ali;
	ali = Malloc( sizeof( ALI ));
	ali->align_i = Malloc ( size * sizeof(int));
	ali->align_j = Malloc ( size * sizeof(int));
	ali->align_s = Malloc ( size * sizeof(CALCTYPE));
	return ali;
}

ARRAYALI *AllocateArrayAliMemory ( size)
int size;
{
	ARRAYALI * ali;
	int i;
	ali = Malloc( sizeof( ARRAYALI ));
	ali->ali    = Malloc ( (size + 1) * sizeof(int));
	ali->scores = Malloc ( (size + 1) * sizeof(CALCTYPE));
	for (i = 0; i <= size; i++) {
		ali->ali[i]   = 0;
		ali->scores[i] = 0;
	}
	ali->length = size;
	return ali;
}

void FreeArrayAliMemory( a )
ARRAYALI *a ;
{
	free( a->ali );
	free( a->scores );
	free( a );
}

COUNTCOLUMN *AllocateCountsMemory ( size )
int size;
{
	COUNTCOLUMN *c = (COUNTCOLUMN*)Malloc( ( size + 1) * sizeof( COUNTCOLUMN ));
	int i, j;
	for (i = 1; i <= size; i++ ) 
		for (j = 0; j < PROFILEWIDTH; j++)
			c[i][j] = 0;

	return (c);
}

FREQUENCYCOLUMN *AllocateFrequenciesMemory ( size )
int size;
{
	FREQUENCYCOLUMN *c = (FREQUENCYCOLUMN*)Malloc( ( size + 1) * sizeof( FREQUENCYCOLUMN ));
	int i, j;
	for (i = 1; i <= size; i++ ) 
		for (j = 0; j < PROFILEWIDTH; j++)
			c[i][j] = 0;

	return (c);
}

PROFILECOLUMN *AllocateProfileMemory ( size )
int size;
{
	PROFILECOLUMN *c = (PROFILECOLUMN*)Malloc( ( size + 1) * sizeof( PROFILECOLUMN ));
	int i, j;
	for (i = 1; i <= size; i++ ) 
		for (j = 0; j < PROFILEWIDTH; j++)
			c[i][j] = 0;

	return (c);
}

void FreeCountsMemory ( c ) 
COUNTCOLUMN *c;
{
	free (c);
}

void FreeProfileMemory ( c ) 
PROFILECOLUMN *c;
{
	free (c);
}

void FreeFrequenciesMemory ( c ) 
FREQUENCYCOLUMN *c;
{
	free (c);
}

void FreeAliMemory ( a ) 
ALI *a;
{
	free(a->align_i);
	free(a->align_j);
	free(a->align_s);
	free (a);
}

void FreeDotsMemory ( dots ) 
DOTS *dots;
{
	free (dots->row);
	free (dots->col);
	free (dots->score);
	free (dots);
}

DOTS* AllocateDotsMemory ( ndots ) 
int ndots;    
{
	DOTS *dots;
	dots = (DOTS*) Malloc( sizeof( DOTS ));
	dots->row   = (int*)Malloc (ndots * sizeof( int )); 
	dots->col   = (int*)Malloc (ndots* sizeof( int ));
	dots->score = (CALCTYPE*)Malloc (ndots * sizeof( CALCTYPE ));
	dots->ndots         = ndots;
	dots->ndotsdiagonal = ndots;

	return (dots);
}

SEQTYPE* CopySequence ( sequence, lsequence ) 
SEQTYPE *sequence;
int lsequence;
{
	SEQTYPE *s = AllocateSequenceMemory( lsequence );
	memcpy( s, sequence, (lsequence + 1) * sizeof (SEQTYPE) );
	return s;
}


COUNTCOLUMN * CopyCounts( counts, lcounts ) 
COUNTCOLUMN * counts;
int lcounts;
{
	COUNTCOLUMN *newcounts;
	newcounts = AllocateCountsMemory( lcounts );

	memcpy( newcounts, counts, (lcounts + 1) * sizeof(COUNTCOLUMN));
	return (newcounts);
}


DOTS *CopyDots( dots ) 
DOTS * dots;
{
	DOTS *newdots;
	newdots = AllocateDotsMemory( dots->ndotsdiagonal );
	memcpy( newdots->row  , dots->row,   dots->ndotsdiagonal * sizeof(int));
	memcpy( newdots->col  , dots->col,   dots->ndotsdiagonal * sizeof(int));
	memcpy( newdots->score, dots->score, dots->ndotsdiagonal * sizeof(CALCTYPE));
	newdots->ndots = dots->ndots;
	newdots->ndotsdiagonal = dots->ndotsdiagonal;

	return (newdots);
}


DOTS *CopyDotsHalf( dots, half ) 
DOTS * dots;
int half;
{
	DOTS *newdots, *cleandots;
	int i;
	int ndots = 0;
	int ndotsdiagonal;

	newdots = AllocateDotsMemory( dots->ndotsdiagonal );

	if (half == UPPER) {
		for (i = 0; i < dots->ndots; i++) {
			if (dots->col[i] > dots->row[i]) {
				newdots->row[ndots]   = dots->row[i];
				newdots->col[ndots]   = dots->col[i];
				newdots->score[ndots] = dots->score[i];
				ndots++;
			}
		}
	} else {
		for (i = 0; i < dots->ndots; i++) {
			if (dots->col[i] < dots->row[i]) {
				newdots->row[ndots]   = dots->row[i];
				newdots->col[ndots]   = dots->col[i];
				newdots->score[ndots] = dots->score[i];
				ndots++;
			}
		}
	}

	ndotsdiagonal = ndots;
	for (i = dots->ndots; i < dots->ndotsdiagonal; i++) {
		newdots->row[ndotsdiagonal]   = dots->row[i];
		newdots->col[ndotsdiagonal]   = dots->col[i];
		newdots->score[ndotsdiagonal] = dots->score[i];
		ndotsdiagonal++;
	}
	newdots->ndots         = ndots;
	newdots->ndotsdiagonal = ndotsdiagonal;

	if (verbose > LL2) 
		printf("Dots reduced to %i (%i)\n", ndots, ndotsdiagonal);


	cleandots = CopyDots( newdots );
	FreeDotsMemory( newdots );
	return (cleandots);
}


/*---------------------------------------------------------------------*/	    
/* Title  : FoldCounts                                                */
/* Funct. : calculates the counts for a profile-alignment       */
/* Author : Andreas Heger		                               */
/* Created: 6.1.1999                                                   */
/*---------------------------------------------------------------------*/	    
COUNTCOLUMN * FoldCounts( origcounts, ali, lrepeat ) 
COUNTCOLUMN *origcounts;
ARRAYALI    *ali;
int          lrepeat;
{
	COUNTCOLUMN   *counts;
	int i, j;

	/* build profile */
	counts = AllocateCountsMemory(  lrepeat );

	for (i = 1; i <= ali->length; i++) {
		if (ali->ali[i] != 0) {
			for (j = 0; j < PROFILEWIDTH; j++) {
				counts[ali->ali[i]][j] += origcounts[i][j];
			}
		}
	}

	return (counts);
}

/*---------------------------------------------------------------------*/	    
/* Title  : FoldDots                                                    */
/* Funct. :                                                             */
/* Author : Andreas Heger		                               */
/* Created: 6.1.1999                                                    */
/*---------------------------------------------------------------------*/	    

/* assumes, that dots are in upper diagonal, d.h. col > row */
/* alignment is to profile, i.e. starts at 1 */

DOTS *FoldDots( dots, origcounts, sequence, lsequence, ali, lrepeat, origprofile ) 
DOTS        *dots;
COUNTCOLUMN *origcounts;
SEQTYPE     *sequence;
int          lsequence;
ARRAYALI    *ali;
int          lrepeat;
PROFILECOLUMN *origprofile;
{
	DOTS *newdots,  *tempdots;
	PROFILECOLUMN   *profile;
	COUNTCOLUMN     *counts;
	FREQUENCYCOLUMN *frequencies, *origfrequencies;

	int i, ndots, row, col;
	int * matrix;

	/* get frequencies and profile */
	profile         = Malloc( (lrepeat + 1)   * sizeof(PROFILECOLUMN));
	frequencies     = Malloc( (lrepeat + 1)   * sizeof(FREQUENCYCOLUMN));
	origfrequencies = Malloc( (lsequence + 1) * sizeof(FREQUENCYCOLUMN));
	counts          = FoldCounts( origcounts, ali, lrepeat);

	Counts2Frequencies( counts,     lrepeat  , frequencies);
	Counts2Frequencies( origcounts, lsequence, origfrequencies);
	Counts2Profile( counts, profile, lrepeat );

	FreeCountsMemory( counts );
	/* build dotplot for sequence-profile-alignment */
	/* row = sequence, col = profile */

	matrix = Malloc( (lrepeat + 1) * (lsequence + 1) * sizeof(int) );
	for (i = 0; i < (lrepeat + 1) * (lsequence + 1); i++) matrix[i] = 0;
	for (i = 0; i < dots->ndots; i++) {
		/* add dot to the alignment, if either one of them is part of the alignment, but not both */
		if ( ((ali->ali[dots->col[i]]) != 0) &&	/* if residue is aligned to a residue in alignment */
				((ali->ali[dots->row[i]] == 0)) ) {	/* and itself is not part of the alignment */
			row = dots->row[i];             /* sequence */
			col = ali->ali[dots->col[i]];   /* profile */
			matrix[ col * lsequence + row] = 1;
			continue;
		} else {
			if ((ali->ali[dots->col[i]] == 0) &&	/* if residue is aligned to a residue in alignment */
					(ali->ali[dots->row[i]] != 0) ) {	/* and itself is not part of the alignment */
				row = dots->col[i];             /* sequence */
				col = ali->ali[dots->row[i]];   /* profile */
				matrix[ col * lsequence + row] = 1;
				continue;
			}
		}
	}


	ndots = 0;
	tempdots = AllocateDotsMemory( dots->ndots);
	for (col = 1; col <= lrepeat; col++) {
		for (row = 1; row <= lsequence; row++) {
			if (matrix[col * lsequence + row] != 0) {
				tempdots->row[ndots] = row;
				tempdots->col[ndots] = col;
				tempdots->score[ndots] = ProfileScore( profile[col],     frequencies[col],         
						origprofile[row], origfrequencies[row]);
				ndots++;
			}
		}
	}

	tempdots->ndots         = ndots;
	tempdots->ndotsdiagonal = ndots;
	free( matrix );

	newdots = CopyDots( tempdots );
	FreeDotsMemory( tempdots );

	FreeProfileMemory( profile );

	free( frequencies );
	free( origfrequencies );

	return (newdots );
}

/*---------------------------------------------------------------------*/	    
/* Title  : CopyDotsWindow                                             */
/* Funct. : copy dots to a new array, but only keep dots inside specified */
/*          window                                                     */
/* Author : Andreas Heger		                               */
/* Created: 11.1.2000                                                  */
/*---------------------------------------------------------------------*/	    

DOTS *CopyDotsWindow( dots, from, to ) 
DOTS * dots;
int from, to;
{
	DOTS *newdots, *tempdots;
	int i, ndots, ndotsdiagonal;
	tempdots = AllocateDotsMemory( dots->ndotsdiagonal );
	ndots = 0;
	for (i = 0; i < dots->ndots; i++) {
		if (dots->row[i] >= from && 
				dots->row[i] <= to   &&
				dots->col[i] >= from &&
				dots->col[i] <= to   ) {
			tempdots->row[ndots]   = dots->row[i];
			tempdots->col[ndots]   = dots->col[i];
			tempdots->score[ndots] = dots->score[i];
			ndots ++;
		}
	}
	ndotsdiagonal = ndots;
	for (i = dots->ndots; i < dots->ndotsdiagonal; i++) {
		if (dots->row[i] >= from && 
				dots->row[i] <= to   &&
				dots->col[i] >= from &&
				dots->col[i] <= to   ) {
			tempdots->row[ndotsdiagonal]   = dots->row[i];
			tempdots->col[ndotsdiagonal]   = dots->col[i];
			tempdots->score[ndotsdiagonal] = dots->score[i];
			ndotsdiagonal ++;
		}
	}

	tempdots->ndots         = ndots;
	tempdots->ndotsdiagonal = ndotsdiagonal;

	newdots = CopyDots( tempdots );
	FreeDotsMemory( tempdots );

	return (newdots);
}

/*---------------------------------------------------------------------*/	    
/* Title  : CreateEnhancedDots2				 */
/* Funct. : uses a profile-alignment to rescore dots and add new dots    */
/*          that are consistent. Both of the dots have to be part of     */
/*          the alignment, if they are to be added. This applies to dots */
/*          within and without the diagonal                              */    
/* Author : Andreas Heger		                                */
/* Created: 6.1.1999                                                     */
/*---------------------------------------------------------------------  */	    
DOTS *CreateEnhancedDots2( origdots, ali, sequence, lsequence, lrepeat, origcounts, origprofile ) 
DOTS * origdots;
ARRAYALI * ali;
SEQTYPE *sequence;
int lsequence;
int lrepeat;
COUNTCOLUMN *origcounts;
PROFILECOLUMN *origprofile;
{
	DOTS * tempdots, * newdots;
	PROFILECOLUMN *profile;
	COUNTCOLUMN * counts;
	FREQUENCYCOLUMN * frequencies, *origfrequencies;
	int ndots, ndotsdiagonal, i, j, from, to, row, col;
	int thisrow, thiscol, lastcol, d, p_row, p_col;
	int *matrix;
	CALCTYPE score;

	int maxdistance = (int)(lrepeat * reducefactor); /* width of bad around diagonal */

	tempdots = CopyDots( origdots );

	/* get frequencies and profile */
	counts      = FoldCounts( origcounts, ali, lrepeat);
	profile     = Malloc( (lrepeat + 1) * sizeof(PROFILECOLUMN));
	frequencies     = Malloc( (lrepeat + 1) * sizeof(FREQUENCYCOLUMN));
	origfrequencies = Malloc( (lsequence + 1) * sizeof(FREQUENCYCOLUMN));

	Counts2Frequencies( counts,     lrepeat  , frequencies);
	Counts2Frequencies( origcounts, lsequence, origfrequencies);
	Counts2Profile( counts, profile, lrepeat );

	/* sort dots by row, do not look at dots on diagonal */
	SortDots( tempdots->row, tempdots->col, tempdots->score, tempdots->ndots);

	/* get aligned region */
	from = 0;         while ((from <= lsequence) && (ali->ali[from] == 0)) from++;
	to   = lsequence; while ((to   > 0         ) && (ali->ali[to]   == 0)) to--;

	if (verbose > LL2) 
		printf("Build enhanced dotarray in range %i-%i with width %i.\n", from, to, maxdistance);

	/* matrix is the band around the diagonal, where dots are added. No dots further than lrepeat apart */
	/* are needed, since we want to reduce the diagonal */
	matrix = Malloc( (lsequence + 1) * maxdistance * sizeof(int)  ); 
	for (i = 0; i < (lsequence + 1)  * maxdistance; i++) matrix[i] = 0;

	/* advance to aligned region */
	i = 0; while (i < tempdots->ndots && tempdots->row[i] < from) i++;

	/* put dots in matrix */
	while (i < tempdots->ndots && tempdots->row[i] <= to) {
		thisrow = tempdots->row[i];

		/* advance to upper diagonal, must exist, since the dots are symmetrical */
		while (thisrow > tempdots->col[i]) i++;

		/* add dots within maxdistance distance to matrix */
		while ( (thisrow == tempdots->row[i]) && 
				(i < tempdots->ndots) &&
				((d = (thiscol = tempdots->col[i]) - thisrow) < maxdistance) ) {

			if ( (ali->ali[thiscol] != 0) && /* only add dots, which are both part of the alignment */ 
					(ali->ali[thisrow] != 0) )
				matrix[thisrow + lsequence * d] = 1;
			i++;
		}

		/* map dots outside this distance to matrix, only if they are both part of the alignment */
		lastcol = tempdots->col[i];
		while ( (thisrow == tempdots->row[i]) && 
				(i < tempdots->ndots) &&
				( (thiscol = tempdots->col[i]) <= to) ) {

			if ( (d = thiscol - lastcol) < maxdistance && 
					ali->ali[thiscol] != 0 &&
					ali->ali[lastcol] != 0
			) 
				matrix[lastcol + d * lsequence] = 1;

			lastcol = thiscol;
			i++;
		}

		while (i < tempdots->ndots && thisrow == tempdots->row[i]) i++;
	}					      

	FreeDotsMemory(tempdots);

	/* build dots, discard dots, which are not part of the alignment */
	tempdots = AllocateDotsMemory( MAXNDOTS );
	ndots = 0;

	for (i = 1; i <= lsequence; i++) {
		for (j = minlrepeat; j < maxdistance; j++) {
			if (matrix[ j * lsequence + i] == 1) {
				row = i; 
				col = i + j;

				p_row = ali->ali[row];
				p_col = ali->ali[col];
				/* alternatively, copy new columns into profile */
				if (p_row != 0 ) {
					if ( p_col != 0) {
						score = ProfileScore( profile[p_row], frequencies[p_row], 
								profile[p_col], frequencies[p_col]);
					} else {
						score = ProfileScore( profile[p_row], frequencies[p_row], 
								origprofile[col], origfrequencies[col] );
					}
				} else {
					if ( p_col != 0) {
						score = ProfileScore( origprofile[row], origfrequencies[row], 
								profile[p_col], frequencies[p_col]);
					} else {
						score = ProfileScore( origprofile[row], origfrequencies[row], 
								origprofile[col], origfrequencies[col] );
					}
				}

				tempdots->row[ndots] = row; tempdots->col[ndots] = col; tempdots->score[ndots] = score;
				ndots ++;
				tempdots->col[ndots] = row; tempdots->row[ndots] = col; tempdots->score[ndots] = score;
				ndots ++;
			}
		}
	}

	free (matrix);
	ndotsdiagonal = ndots;

	/* add diagonal; if dot was part of the alignment, use the profile-score, otherwise copy score from old alignemt */
	i = origdots->ndots; while (origdots->row[i] < from) i++;

	while (origdots->row[i] <= to && i < origdots->ndotsdiagonal) {
		thisrow = origdots->row[i];
		tempdots->row[ndotsdiagonal] = thisrow;
		tempdots->col[ndotsdiagonal] = thisrow;
		p_row = ali->ali[thisrow];
		if (p_row != 0) 
			tempdots->score[ndotsdiagonal] = ProfileScore( profile[p_row], frequencies[p_row], profile[p_row], frequencies[p_row]);
		else
			tempdots->score[ndotsdiagonal] = ProfileScore( origprofile[thisrow], origfrequencies[thisrow], 
					origprofile[thisrow], origfrequencies[thisrow] );
		ndotsdiagonal++;
		i++;
	}

	FreeProfileMemory( profile );
	FreeCountsMemory ( counts );

	free( frequencies );
	free( origfrequencies );

	tempdots->ndots = ndots;
	tempdots->ndotsdiagonal = ndotsdiagonal;

	newdots = CopyDots( tempdots );
	FreeDotsMemory( tempdots );

	return (newdots );

}



ALI *CopyAli( ali ) 
ALI * ali;
{
	ALI *newali;
	int newsize = ali->lastindex + 1;
	newali = AllocateAliMemory( newsize);
	memcpy( newali->align_i, ali->align_i, newsize * sizeof(int));
	memcpy( newali->align_j, ali->align_j, newsize * sizeof(int));
	memcpy( newali->align_s, ali->align_s, newsize * sizeof(CALCTYPE));
	newali->length    = ali->length;
	newali->lastindex = ali->lastindex;
	newali->score     = ali->score;

	return (newali);
}

ARRAYALI *CopyArrayAli( ali )
ARRAYALI *ali;
{
	ARRAYALI *newali;
	int newsize = ali->length + 1;
	newali = AllocateArrayAliMemory( ali->length );
	memcpy( newali->ali,    ali->ali,    newsize * sizeof(int));
	memcpy( newali->scores, ali->scores, newsize * sizeof(CALCTYPE));
	newali->length = ali->length;

	return (newali);
}


DOTS* MaskDots ( dots, 
		sequence, 
		lsequence
)
DOTS *dots;
SEQTYPE  *sequence;
int      lsequence;
{ 

	int i, x;
	DOTS *tempdots, *newdots;

	tempdots = AllocateDotsMemory( dots->ndotsdiagonal ); 

	if (verbose > LL3) {
		printf("MaskDots called with %i: \n", lsequence);
		printf("Old number of dots: %i %i\n", dots->ndots, dots->ndotsdiagonal);
		PrintDots( dots->row, dots->col, dots->score, dots->ndotsdiagonal );
		PrintSequence( sequence, lsequence);
	}

	x = 0;
	/* copy non-diagonal */
	for (i = 0; i < dots->ndots; i++) {
		if ( (sequence[dots->row[i]] != MASKCODE) &&
				(sequence[dots->col[i]] != MASKCODE) ) {
			tempdots->row[x]   = dots->row[i];
			tempdots->col[x]   = dots->col[i];
			tempdots->score[x] = dots->score[i];
			x++;
		}
	}
	tempdots->ndots = x;
	/* copy diagonal */
	for (i = dots->ndots; i < dots->ndotsdiagonal; i++) {
		if ( (sequence[dots->row[i]] != MASKCODE) &&
				(sequence[dots->col[i]] != MASKCODE) ) {
			tempdots->row[x]   = dots->row[i];
			tempdots->col[x]   = dots->col[i];
			tempdots->score[x] = dots->score[i];
			x++;
		}
	}
	tempdots->ndotsdiagonal = x;
	newdots = CopyDots( tempdots );

	FreeDotsMemory( tempdots );

	if (verbose > LL3) {
		printf("New number of dots: %i %i\n", newdots->ndots, newdots->ndotsdiagonal);
		PrintDots( newdots->row, newdots->col, newdots->score, newdots->ndotsdiagonal);
	}

	return (newdots);
}


ARRAYALI * Ali2ArrayAli ( ali, lsequence ) 
ALI * ali;
{
	int i;
	ARRAYALI *newali = AllocateArrayAliMemory( lsequence );

	/* changed, since I changed TraceBackDotAlignment */
	for (i = 0; i <= ali->lastindex; i++) {
		newali->ali[ali->align_j[i]]    = ali->align_i[i];
		newali->scores[ali->align_j[i]] = ali->align_s[i];
	}
	newali->length = lsequence;
	newali->score  = ali->score;

	return (newali);
}

/* Find best window in tubealignment */
/* assumes, that the alignment is backwards, i.e. 
   1. starting at the high residue numbers 
   2. align_i > align_j
   3. first repeat is unaligned

   for example | 6 4 x| 5 3 x | 4 2 x | 3 1 x |
 */

/* assumes that the alignment is forwards, i.e. 
   1 2 3 4 5  6 7 8 9 10
       7 8 9 10
 */
/*--------------------------------------------------------------------- */	    
/* Title  : FindBestWindowPairWise                                      */
/* Funct. : search for best window in alignment with length of maxdiag  */
/*          (note: the window-length is on residues, not on alignment)  */
/* Author : Andreas Heger		                               */
/* Created: 6.1.1999                                                    */
/*--------------------------------------------------------------------- */	    
SEQUENCEWINDOW FindBestWindowPairWise ( sequence, lsequence, aali, maxdiag)
SEQTYPE *sequence;
int lsequence;
ARRAYALI *aali;
int maxdiag;
{
	int i,j,x, lasti, lastj, besti, bestj;
	CALCTYPE s  = 0;
	CALCTYPE bests = 0;
	SEQUENCEWINDOW bestwindow; 

	i = 1; while (i < aali->length && aali->ali[i] == 0) i++; /* i = trailing index on sequence */
	j = i + maxdiag - 1;		                      /* j = leading  index on sequence */
	x = aali->length; while (x > 0 && aali->ali[x] == 0) x--; /* x = end of aligned region */

	besti = i;
	bestj = j;

	lasti = i++;
	lastj = j++;
	while ( j <= x ) { /* remember we are counting down */

		/* case distinction for i*/
		if ( aali->ali[lasti] != 0) {             /* if last was a match */
			s -= aali->scores[lasti];              /* substract score */
		} else {
			if ( aali->ali[i] == 0 )  	          /* if it was a gap */ 
				s -= bw_gapelon_i;                /* reduce gap penalty for elongation */
			else 
				s -= (bw_gapopen_i + bw_gapelon_i); /* reduce gap penalty for opening */
		}

		/* case distinction for j*/
		if ( aali->ali[j] != 0) {             /* if new is a match */
			s += aali->scores[j];              /* add score */
		} else {
			if ( aali->ali[lastj] == 0 )  	  /* if it was a gap */ 
				s += bw_gapelon_i;                /* add gap penalty for elongation */
			else 
				s += (bw_gapopen_i + bw_gapelon_i); /* add gap penalty for opening */
		}

		if (s > bests) {
			besti = i; bestj = j; bests = s;
		}

		/* printf("%5i %5i %5i %5i %5i %5i %5i %5i %5.2f %5.2f\n", i, j, aali->ali[i], aali->ali[j],  
	       lasti, lastj, aali->ali[i], aali->ali[j], 
	       aali->scores[i], aali->scores[j]); */
		lasti = i++;
		lastj = j++;
	}

	if (bestj > lsequence) /* bug-fix, in future calculate best window more elegantly with variable length */
		bestj = lsequence;

	bestwindow.from = besti;
	bestwindow.to   = bestj;

	bestwindow.score= s;

	return (bestwindow);
}    
/*--------------------------------------------------------------------- */    	    
/* Title  : PruneWindow                                                    */
/* Funct. : remove ends from a window, that contribute to a negative score */
/* Author : Andreas Heger		                               */
/* Created: 15.1.1999                                                    */
/*--------------------------------------------------------------------- */	    
SEQUENCEWINDOW PruneWindow( window, ali ) 
SEQUENCEWINDOW window;
ARRAYALI *ali;
{
	int i;
	SEQUENCEWINDOW newwindow;
	i = window.from;
	while ( (i <= ali->length) && (ali->ali[i] == 0 || ali->scores[i] < 0) ) i++;
	if (i <= ali->length) 
		newwindow.from = i;
	else 
		newwindow.from = 0;

	i = window.to;
	while ( (i > 0)  && (ali->ali[i] == 0 || ali->scores[i] < 0) ) i--;
	if (i >= newwindow.from) 
		newwindow.to = i;
	else
		newwindow.to = 0;

	return (newwindow);
}


/*--------------------------------------------------------------------- */    	    
/* Title  : GetBestWindow                                                  */
/* Funct. : remove ends from a window, that contribute to a negative score */
/* Author : Andreas Heger		                               */
/* Created: 15.1.1999                                                    */
/*--------------------------------------------------------------------- */	    
SEQUENCEWINDOW GetBestWindow ( from, to, maxdiag, sequence, lsequence, dots)
int from, to;
int maxdiag;
SEQTYPE *sequence;
int lsequence;
DOTS *dots;

{
	int d1;
	ALI *ali;
	ARRAYALI *aali;
	SEQUENCEWINDOW bestwindow;
	DOTS *tempdots = CopyDots( dots ); /* copying is necessary, since alignment routine sorts dots */

	bestwindow.from = 0;
	bestwindow.to = 0;

	d1 = rint(maxdiag / scalefactor); 

	if (d1 < minsamplediag) d1 = minsamplediag; 

	ali = AllocateAliMemory( MAXALI_SIZE);
	dotalign_tube(   tempdots->row, tempdots->col, tempdots->score, tempdots->ndots,
			0,							/* no roll */ 
			from, to, 
			from, to, 
			rd_gapopen_i, rd_gapelon_i,
			rd_gapopen_j, rd_gapelon_j,
			maxdiag - d1, 
			maxdiag + d1,
			ali );
	FreeDotsMemory( tempdots );

	if (verbose > LL1) 
		printf ("Score of tubealignment: %5.2f in tube (%i;%i)\n",
				ali->score, maxdiag - d1, maxdiag + d1);

	if (verbose > LL2)
		DumpAlignment(ali, "FindBestWindow");

	if ( ali->score < mintubescore) {
		FreeAliMemory( ali );
		return bestwindow;
	}

	aali = Ali2ArrayAli( ali, lsequence); 
	FreeAliMemory( ali );
	bestwindow = FindBestWindowPairWise( sequence, lsequence, aali, maxdiag); 

	if (verbose > LL1) 
		printf ("BestWindow before pruning on diagonal %i: %i to %i\n", maxdiag, bestwindow.from, bestwindow.to);

	bestwindow = PruneWindow( bestwindow, aali );
	if (verbose > LL1) 
		printf ("BestWindow after pruning on diagonal  %i: %i to %i\n", maxdiag, bestwindow.from, bestwindow.to);

	FreeArrayAliMemory( aali );

	return (bestwindow);
}

/*---------------------------------------------------------------------*/	    
/* Title  : GetDiagonal                                                */
/* Funct. : retrieves the highest scoring diagonal                     */
/*          that is further than minlrepeat from the main diagonal apart */
/* Author : Andreas Heger		                               */
/* Created: 4.1.1999                                                   */
/*---------------------------------------------------------------------*/	    

int GetDiagonal ( from, to, dots, bestscore, iteration_start ) 
int from, to;
DOTS *dots;
CALCTYPE *bestscore;
int iteration_start;
{
	int maxdiag = 0;
	int iteration = iteration_start;

	while (  ((maxdiag    < minlrepeat) || (*bestscore < mindiagscore) ) 
			&& iteration < MAX_NITERATIONS ) { 
		maxdiag = getnbestdiagpostrace( dots->row, dots->col, dots->score, dots->ndots,
				from, to, from, to, 
				iteration, /* get best level */
				diagwidth,
				-RINFINITE, /* use positive and negative scores */
				bestscore );


		if (maxdiag < 0) maxdiag = -maxdiag; 

		if (verbose > LL3)
			printf ("Maximum Diagonal (score): %i (%5.2f)\n", maxdiag, *bestscore);

		iteration ++;

		if (maxdiag == 0 || *bestscore == 0)
			break;
	}

	/* if we have not found a diagonal, try to find the highest scoring diagonal taking only into account positive scores */
	/* this is useful, if you have only short repeats, that are obscured by low/negative scoring traces in the same diagonal */

	if (iteration == MAX_NITERATIONS) {
		if (verbose > LL1) 
			printf ("Repeated diagonal calculation with only positive scores\n"); 
		iteration = 0;
		maxdiag   = 0;

		while ( maxdiag < minlrepeat && iteration < MAX_NITERATIONS) { 
			maxdiag = getnbestdiagsum( dots->row, dots->col, dots->score, dots->ndots,
					from, to, from, to, 
					iteration, /* get best level */
					diagwidth,
					0.0, /* use positive scores*/
					bestscore );


			if (maxdiag < 0) maxdiag = -maxdiag; 

			if (verbose > LL3)
				printf ("Maximum Diagonal (score): %i (%5.2f)\n", maxdiag, *bestscore);

			iteration ++;

			if (maxdiag == 0 || *bestscore == 0)
				break;
		}

	}

	if (iteration == MAX_NITERATIONS) 
		maxdiag = 0;

	if (verbose > LL1)
		printf ("Maximum Diagonal (score, iteration): %i (%5.2f, %i)\n", maxdiag, *bestscore, iteration - 1);

	return (maxdiag);
}

/*SEQUENCEWINDOW **/
void PrintPrettyMA ( ali, nrepeats, repeatscores, zscores, ssequence, lsequence, buffer ) 
ARRAYALI *ali;
int nrepeats;
CALCTYPE *repeatscores;
CALCTYPE *zscores;
SEQTYPE * ssequence;
int lsequence;
char * buffer;
{
	int i, j;
	int maxcol = 0;
	char * ma;
	SEQUENCEWINDOW *window;
	int nins, totalinserts;
	int *position;
	int start;
	int nrep;
	int ma_len;
	int insertcol, thiscol, lastcol;
	int first, lasti;
	char *sequence = Sequence2String(ssequence, lsequence); /* note: numbering in string start from 0!!! */

	int *inserts = Malloc( sizeof(int) * lsequence);
	for (i = 0; i < lsequence; i++) inserts[i] = 0;

	i = 1; while ((ali->ali[i] < 1) && i <= (ali->length)) i++;

	start = i;
	lastcol = ali->ali[i];
	/* get maximum number of inserts after each profile position */
	while (i <= ali->length) {
		nins = 0;
		while ((i <= ali->length) && (ali->ali[i] == 0)) { i++; nins++; };
		if ( i > ali->length) break;
		if ( ( ali->ali[i] > lastcol ) && (inserts[ali->ali[i] - 1] < nins )) /* only save for repeat-internal gaps */
			/* 	    inserts[lastcol] = nins; */
			inserts[ali->ali[i] - 1] = nins; /* give gaps to previos profilecolumn, since there might be several profilecolumns missing*/
		lastcol = ali->ali[i];
		if (lastcol > maxcol) maxcol = lastcol;
		i++;
	}

	totalinserts = 0; for (i = 1; i < maxcol; i++) totalinserts += inserts[i];

	position  = (int*) Malloc( sizeof(int) * (maxcol + 1));
	ma_len    = (maxcol + totalinserts + 1); /* +1 for \0 at the end */
	ma        = (char*)Malloc( sizeof(char) * ma_len * nrepeats);
	window    = (SEQUENCEWINDOW*) Malloc( sizeof(SEQUENCEWINDOW) * nrepeats);

	for (i = 0; i < (ma_len * nrepeats); i++) ma[i] = GAPCHAR;
	for (i = ma_len - 1; i < (ma_len * nrepeats); i+= ma_len) ma[i] = '\0'; 

	position[0] = -1;                        /* so that position[1] = 0 */
	inserts[0]  = 0;                         /* calculate positions for each column */
	for (i = 1; i <= maxcol; i++) position[i] = 1 + position[i-1] + inserts[i-1];

	i     = start; 
	lasti = start;
	nrep = 0;                 /* number of repeat */

	lastcol = 0;
	first = 1;
	window[nrep].from = start;
	while (i <= ali->length ) {

		thiscol = ali->ali[i];

		if ( (lastcol > thiscol) || 
				((i - lasti - 1) > splitdistance) ) { /* start a new repeat, if wrapped around profile or separated by large gap */
			window[nrep].to = lasti;
			nrep++; /* new repeat started */
			window[nrep].from = i;
		}

		insertcol = position[thiscol] + ma_len * nrep; 
		/* printf("%i %i %i %i %i %i \n", i, insertcol, thiscol, nrep, ma_len, position[thiscol]); */
		ma[insertcol] = sequence[i-1]; /* insert match */
		if (!first) {
			if (lastcol < thiscol) {          /* if internal gap */
				j = i - 1;
				while (ali->ali[j] == 0) {             /* go backwards for gaps */
					ma[--insertcol] = sequence[--j] + 32; /* convert to lower case */
				}
			} 
		} else {
			first = 0;
		}
		lasti   = i;
		i++;
		lastcol = thiscol;
		while (i <= ali->length && ali->ali[i] == 0) i++; /* jump over gap regions */
	}
	window[nrep].to = lasti;
	for (i = 0; i < nrepeats; i++) {
		sprintf(&buffer[strlen(buffer)],
				"%5i-%5i (%5.2f/%5.2f)\t%s\n", window[i].from, window[i].to, repeatscores[i], zscores[i], &ma[ma_len * i]); 
	}
	free(position);
	free(inserts);
	free(ma);
	free(sequence);

	free(window);

	/* return(window);*/
}

/* Print a pretty alignment in basic form */ 
void PrintPrettyAli2 ( ali, ssequence, lsequence ) 
ARRAYALI *ali;
SEQTYPE * ssequence;
int lsequence;
{
	int i, j;
	int maxcol = 0;
	char * ma;
	SEQUENCEWINDOW *window;
	int nins, totalinserts;
	int *position;
	int start;
	int nrep;
	int ma_len;
	int insertcol, thiscol, lastcol;
	int first, lasti;
	char *sequence = Sequence2String(ssequence, lsequence); /* note: numbering in string start from 0!!! */

	int *inserts = Malloc( sizeof(int) * lsequence);
	for (i = 0; i < lsequence; i++) inserts[i] = 0;

	i = 1; while ((ali->ali[i] < 1) && i <= (ali->length)) i++;

	start = i;
	lastcol = ali->ali[i];
	/* get maximum number of inserts after each profile position */
	while (i <= ali->length) {

		nins = 0;
		while ((i <= ali->length) && (ali->ali[i] == 0)) { i++; nins++; };
		if ( i > ali->length) break;
		if ( ( ali->ali[i] > lastcol ) && (inserts[ali->ali[i] - 1] < nins )) /* only save for repeat-internal gaps */
			/* 	    inserts[lastcol] = nins; */
			inserts[ali->ali[i] - 1] = nins; /* give gaps to previos profilecolumn, since there might be several profilecolumns missing*/
		lastcol = ali->ali[i];
		if (lastcol > maxcol) maxcol = lastcol;
		i++;
	}

	totalinserts = 0; for (i = 1; i < maxcol; i++) totalinserts += inserts[i];

	position  = (int*) Malloc( sizeof(int) * (maxcol + 1));
	ma_len    = (maxcol + totalinserts + 1); /* +1 for \0 at the end */
	ma        = (char*)Malloc( sizeof(char) * ma_len * MAXNREPEATS);
	window    = (SEQUENCEWINDOW*) Malloc( sizeof(SEQUENCEWINDOW) * MAXNREPEATS);

	for (i = 0; i < (ma_len * MAXNREPEATS); i++) ma[i] = GAPCHAR;
	for (i = ma_len - 1; i < (ma_len * MAXNREPEATS); i+= ma_len) ma[i] = '\0'; 

	position[0] = -1;                        /* so that position[1] = 0 */
	inserts[0]  = 0;                         /* calculate positions for each column */
	for (i = 1; i <= maxcol; i++) position[i] = 1 + position[i-1] + inserts[i-1];

	i     = start; 
	lasti = start;
	nrep = 0;                 /* number of repeat */

	lastcol = 0;
	first = 1;
	window[nrep].from = start;
	while (i <= ali->length ) {

		thiscol = ali->ali[i];

		if ( (lastcol > thiscol) || 
				((i - lasti - 1) > splitdistance) ) { /* start a new repeat, if wrapped around profile or separated by large gap */
			window[nrep].to = lasti;
			nrep++; /* new repeat started */
			window[nrep].from = i;
		}

		insertcol = position[thiscol] + ma_len * nrep; 
		/* printf("%i %i %i %i %i %i \n", i, insertcol, thiscol, nrep, ma_len, position[thiscol]); */
		ma[insertcol] = sequence[i-1]; /* insert match */
		if (!first) {
			if (lastcol < thiscol) {          /* if internal gap */
				j = i - 1;
				while (ali->ali[j] == 0) {             /* go backwards for gaps */
					ma[--insertcol] = sequence[--j] + 32; /* convert to lower case */
				}
			} 
		} else {
			first = 0;
		}
		lasti   = i;
		i++;
		lastcol = thiscol;
		while (i <= ali->length && ali->ali[i] == 0) i++; /* jump over gap regions */
	}
	window[nrep].to = lasti;

	for (i = 0; i <= nrep; i++) {
		printf("%5i-%5i \t%s\n", window[i].from, window[i].to, &ma[ma_len * i]); 
	}
	free(position);
	free(inserts);
	free(ma);
	free(sequence);
	free(window);

	return;
}



/*--------------------------------------------------------------------- */	    
/* Title  : ReduceProfileAlignment                                      */
/* Funct. : renumbers profile-columns, so that there no empty ones      */
/* Author : Andreas Heger		                               */
/* Created: 9.1.1999                                                    */
/*--------------------------------------------------------------------- */	    

int ReduceProfileAlignment( ali, lrepeat, mincounts )
ARRAYALI * ali;
int lrepeat;
int mincounts; /* minimum number of counts required for keeping a profile-column */
{
	int returnlrepeat;
	int * colcounts = Malloc( sizeof(int) * (lrepeat + 1));
	int * map       = Malloc( sizeof(int) * (lrepeat + 1));
	int i, j;
	int mincutoff;
	int maxcutoff;

	for (i = 0; i <= lrepeat; i++) colcounts[i] =0;
	for (i = 1; i <= ali->length; i++) colcounts[ali->ali[i]] ++;

	/* get core region of alignment */
	mincutoff = 1;       while (colcounts[mincutoff] < mincounts) mincutoff++;
	maxcutoff = lrepeat; while (colcounts[maxcutoff] < mincounts) maxcutoff--;

	/* create a map for skipping of empty columns */
	i = 1; j = 0;
	map[0] = 0;
	while (i <= lrepeat ) {
		if (colcounts[i] >= mincounts && i >= mincutoff && i <= maxcutoff ) 
			map[i] = ++j;
		else
			map[i] = 0;
		i++;
	}

	/* save the new accurate repeat length */
	returnlrepeat = j;

	/* iterate through alignment and apply map */
	for (i = 1; i <= ali->length; i++) ali->ali[i] = map[ali->ali[i]];

	free( map );
	free( colcounts );

	return (returnlrepeat);
}

/* Score Alignment against itself */
/* offset: for resetting alignment to 1 */
/* use 0, if you have a profile-alignment */
void AutoScore( ali, lrepeat, offset, sequence )
ARRAYALI * ali;
int lrepeat;
int offset;
SEQTYPE * sequence;
{

	COUNTCOLUMN    *counts  = AllocateCountsMemory (lrepeat); 
	PROFILECOLUMN  *profile = AllocateProfileMemory(lrepeat);
	int i;

	/*     DumpArrayAlignment( ali, "Autoscore"); */
	/* calculate a profile and score repeats against profile */
	for (i = 1; i <= ali->length; i++) 
		if (ali->ali[i] != 0) 
			counts[ali->ali[i] - offset][sequence[i]]++;

	Counts2Profile( counts, profile, lrepeat);

	for (i = 1; i<= ali->length; i++) {
		if (ali->ali[i] != 0)
			ali->scores[i] = profile[ali->ali[i] - offset][sequence[i]];
	}

	FreeCountsMemory( counts );
	FreeProfileMemory( profile );
}


/*----------------------------------------------------------------------------------------*/
ARRAYALI * SelectRepeatsStatistical ( ali, sequence, lsequence, lrepeat, gapopen, gapelon,
		repeatscores, zscores, returnnrepeats, returnlrepeat
)
ARRAYALI * ali;
SEQTYPE *sequence;
int lsequence, lrepeat;
CALCTYPE gapopen, gapelon;
CALCTYPE *repeatscores;
CALCTYPE *zscores;
int *returnnrepeats;
int *returnlrepeat;

{

	SEQUENCEWINDOW *windows = Malloc( MAXNREPEATS * sizeof(SEQUENCEWINDOW));
	COUNTCOLUMN    *counts  = AllocateCountsMemory(lrepeat);
	PROFILECOLUMN  *profile = AllocateProfileMemory(lrepeat);
	ARRAYALI *newali = AllocateArrayAliMemory( ali->length );
	int nrepeats, i, j;
	int nsavedrepeats = 0;
	double mean, std, zscore;

	newali->length = ali->length;

	for (i = 0; i <= ali->length; i++)
		newali->ali[i] = 0;

	/* calculate a profile and score repeats against profile------------------------------------------------------------- */
	/* build counts */
	for (i = 1; i <= ali->length; i++) 
		if (ali->ali[i] != 0) 
			counts[ali->ali[i]][sequence[i]]++;

	/* calculate profile */
	Counts2Profile( counts, profile, lrepeat);

	/* rescore */
	for (i = 1; i<= ali->length; i++) {
		if (ali->ali[i] != 0)
			ali->scores[i] = profile[ali->ali[i]][sequence[i]];
	}

	/* count and score repeats */
	nrepeats = CountAndScoreRepeats( ali, lrepeat, gapopen, gapelon, windows );

	/* calculate z-scores for each repeat ------------------------------------------------------------------------- */
	/* get mean and std-deviation for shuffled sequence */
	CalculateZScore( sequence, lsequence,
			profile, lrepeat,
			gapopen, gapelon,
			gapopen, gapelon,
			shuffle_iterations,
			&mean, &std );
	if (verbose > LL2) 
		printf("SelectRepeatsStatistical: %5.2f +/- %5.2f for %i iterations\nZScores: ", mean, std, shuffle_iterations) ;

	for (i = 0; i < nrepeats; i++) {
		zscore = (windows[i].score - mean) / std;
		if (verbose > LL2) printf("%5.2f ", zscore);
		if (zscore > minzscore) {
			for (j = windows[i].from; j <= windows[i].to; j++) 
				newali->ali[j] = ali->ali[j];
			zscores[nsavedrepeats] = zscore;
			nsavedrepeats++;
		}
	}

	if (verbose > LL2) printf("\n");

	/* set correct scores in new alignment */
	for (i = 1; i<= newali->length; i++) {
		if (newali->ali[i] != 0)
			newali->scores[i] = profile[newali->ali[i]][sequence[i]];
		else
			newali->scores[i] = 0;
	}

	/* cut away overhanging ends and close gaps in MA, only sensible for MA with more than one sequence*/
	if (nsavedrepeats > 1)			/* note: I have to return a single repeat for masking */
		lrepeat = ReduceProfileAlignment( newali, lrepeat, 2 );

	/* get the final scores for the repeats */
	nrepeats = CountAndScoreRepeats( newali, lrepeat, gapopen, gapelon, windows );
	for (i = 0; i < nrepeats; i ++) 
		repeatscores[i] = windows[i].score;

	*returnnrepeats = nrepeats;
	*returnlrepeat  = lrepeat;

	free( windows );
	FreeCountsMemory( counts );
	FreeProfileMemory( profile );

	return (newali);
}  

/*--------------------------------------------------------------------------*/
void WellcomeMsg() 
{

	if (verbose > LL0) {
		/* printf("------ Supplied Parameters ------------------------------------------------\n");
		printf("%s %s %s\n", filename_ma, filename_sequence, filename_lfasta); */
		printf("------ Internal Parameters ------------------------------------------------\n");
		printf("mindiag      : %5i\tmaxlevels    : %5i\tdiagwidth    : %5i\n", mindiag, maxlevels, diagwidth);
		printf("scalefactor  : %5i\tminsamplediag: %5i\tminlrepeat   : %5i\n", scalefactor, minsamplediag, minlrepeat);
		printf("minndots     : %5i\tselectfactor : %5.2f\tmintubescore : %5.2f\n", minndots, selectfactor, mintubescore);
		printf("splitdistance: %5i\tshuffle_iterations: %i\tmaxdiaognals : %i\n", splitdistance, shuffle_iterations, maxdiagonals);
		printf("minswscore   : %5.2f\tcutoffscore  : %5.2f\tloglevel     : %5i\n", minswscore, cutoffscore, verbose);
		printf("mindiagscore : %5.2f\treducefactor : %5.2f\tmincollectscore: %5.2f\n", mindiagscore, reducefactor, mincollectscore);     
		printf("minzscore    : %5.2f\n", minzscore);

		printf("Gap-Penalties (Diagonal Selection   ): %6.2f %6.2f %6.2f %6.2f\n", rd_gapopen_i, rd_gapelon_i, rd_gapopen_j, rd_gapelon_j);
		printf("Gap-Penalties (Bestwindow           ): %6.2f %6.2f %6.2f %6.2f\n", bw_gapopen_i, bw_gapelon_i, bw_gapopen_j, bw_gapelon_j);
		printf("Gap-Penalties (Dotalignment         ): %6.2f %6.2f %6.2f %6.2f\n", gapopen_i, gapelon_i, gapopen_j, gapelon_j);
		printf("Gap-Penalties (Collecting repeats   ): %6.2f %6.2f %6.2f %6.2f\n", hp_gapopen_i, hp_gapelon_i, hp_gapopen_j, hp_gapelon_j);
		printf("Gap-Penalties (Register             ): %6.2f %6.2f\n", rg_gapopen, rg_gapelon);
		printf("---------------------------------------------------------------------------\n");
	}

}

/*----------------------------------------------------------------------------------------------*/
/* Title: AlignBestWindow									*/
/* Function: Wrapper for dotalign								*/
/*
   15.3.2001:	Added a copy step for the dots in origdots, so that they are not changed in any way
		during this routine. During the alignmentstep, the dots get sorted. There will be hard 
		to trace bugs, because ndots does not point to the end of the dots without the diagonal.

 */
/*----------------------------------------------------------------------------------------------*/

ARRAYALI * AlignBestWindow( window, lsequence, origdots ) 
SEQUENCEWINDOW window;
int lsequence;
DOTS * origdots;

{
	ALI * ali;
	ARRAYALI * aali;

	DOTS * tempdots = CopyDots( origdots );

	ali      = AllocateAliMemory( MAXALI_SIZE);

	dotalign(   tempdots->row, tempdots->col, tempdots->score, tempdots->ndotsdiagonal,
			1,
			1, lsequence, 
			window.from, window.to, 
			gapopen_i, gapelon_i,
			gapopen_j, gapelon_j, 
			ali );

	aali = Ali2ArrayAli( ali, lsequence );
	FreeAliMemory( ali );
	FreeDotsMemory( tempdots );

	return (aali);
}

/*----------------------------------------------------------
 Title: MapAlignment
 Function: Shifts an alignment, so that residue startcol is
 in the first column with index start.
----------------------------------------------------------*/

void MapAlignment ( startcol, 
		mincol, maxcol, 
		ali,
		start)
int startcol, mincol, maxcol;
ARRAYALI *ali;
int start;
{
	int i, j;
	int lali = ali->length;
	int * map = (int*)Malloc((maxcol + 1) * sizeof(int));

	/* build map */
	j = start - 1;
	for (i = 0; i <= maxcol; i++) map[i] = 0;
	for (i = startcol; i <= maxcol;   i++) map[i]=++j; 
	for (i = mincol;   i <  startcol; i++) map[i]=++j; 

	/* map */
	for (i = 1; i <= lali; i++) ali->ali[i] = map[ali->ali[i]];

	free( map );
	return;
}


/*---------------------------------------------------------------------*/	    
/* Title  : MaskAlignment                                              */
/* Funct. : mask out dots from residues, that are aligned in an        */ 
/*          alignment                                                  */
/* Author : Andreas Heger		                               */
/* Created: 4.1.1999                                                   */
/*---------------------------------------------------------------------*/	    
DOTS * MaskAlignment ( dots, ali ) 
DOTS * dots;
ARRAYALI * ali;
{

	SEQTYPE * mask = AllocateSequenceMemory( ali->length );
	int i;

	if (verbose > LL3) 
		printf("MaskAlignment called\n");

	for (i = 1; i <= ali->length; i++) 
		if (ali->ali[i] != 0) 
			mask[i] = MASKCODE;
		else
			mask[i] = 0;

	mask[0] = MASKCODE;

	return MaskDots( dots, mask, ali->length);
}

/*---------------------------------------------------------------------*/	    
/* Title  : CalculateRepeatMA                                          */
/* Funct. : align repeats with wrapping in one region and then proceed */ 
/*          to the next region by aligning bestwindow to the sequence  */
/* Author : Andreas Heger		                               */
/* Version: 1.1 */
/* Created: 4.1.1999                                                   */
/*---------------------------------------------------------------------*/	    
ARRAYALI * CalculateRepeatMA (    from1, to1, 
		sequence, lsequence,
		origprofile, lprofile, 
		origdots, origcounts,
		level,
		returnscore,
		returnnrepeats,
		returnalignment,
		reducedots)
		int from1, to1;
SEQTYPE *sequence;
int lsequence;
PROFILECOLUMN *origprofile;
int lprofile;
DOTS * origdots;
COUNTCOLUMN * origcounts;
int level;
CALCTYPE *returnscore;
int *returnnrepeats;
char *returnalignment;
DOTS *reducedots;
{

	DOTS *tempdots, *tempdots2;
	SEQUENCEWINDOW bestwindow;
	int lrepeat, lrepeat_new;
	ALI *ali;
	ARRAYALI *intermediate_ali, *smaller_ali, *aali;
	CALCTYPE *repeatscores, *zscores;
	int nrepeats, i, j;
	CALCTYPE sumscore, smaller_score, diagscore;
	int bestdiagonal, smaller_nrepeats;
	int from2, to2,  niterations;
	char buffer[100];
	char resultbuffer[MAXRESULTSIZE];
	char * smaller_alignment;
	int level_diagonal;

	if (verbose > LL1) 
		printf("Entering recursive repeat calculation on level %i (%i/%i) with hopping in window %i-%i\n", 
				level + 1, lsequence, reducedots->ndots, from1, to1);

	*returnscore    = 0;
	*returnnrepeats = 0;

	if (lsequence < (minlrepeat * 2)) return NULL;

	if (reducedots->ndots < minndots) {
		if( verbose > LL1) 
			printf("--->Leaving, because not enough dots\n"); 
		return NULL;
	}

	if (verbose > LL4)
		PrintDots( reducedots->row, reducedots->col, reducedots->score, reducedots->ndotsdiagonal);

	if (verbose > LL3) { 
		sprintf(buffer, "reducedots_%i_%i.png", iteration, level + 1);
	}

	intermediate_ali = smaller_ali = aali = NULL;

	/* Iteration for finding best diagonal and best window ----------------------------------------- */
	level_diagonal = 0;
	bestdiagonal   = 0;
	bestwindow.from= 0;
	lrepeat        = 0;
	diagscore      = 1;   
	while ( level_diagonal < maxdiagonals &&		/* iterate for a maximum number of iterations */
			diagscore != 0 && (				/* and while the the diagonal score is not zero and while */
					bestdiagonal   < minlrepeat   ||		/* we have not found a diagonal with a minimum distance */
					lrepeat        < minlrepeat   ||		/* or we have not a repeat with a minimum length           */
					diagscore      < mindiagscore ||		/* or the minimum diagonal score is below our fancy */
					bestwindow.from == 0) ) {			/* or we couldn't find a best window */

		/* get highest scoring diagonal --------------------------------------------- */

		bestdiagonal = GetDiagonal( 1, lsequence, reducedots, &diagscore, level_diagonal ); 

		if ( bestdiagonal < minlrepeat ) {
			if (verbose > LL1) 
				printf("-----> Iterating, because diagonal too close to main diagonal (%i/%i)\n", bestdiagonal, minlrepeat);
			level_diagonal++;
			continue;
		}

		if (diagscore < mindiagscore ) {
			if (verbose > LL1) 
				printf("-----> Iterating, because diagonal score is too low (%5.2f/%5.2f)\n", diagscore, mindiagscore);
			level_diagonal++;
			continue;
		}

		/* get a template window ---------------------------------------------------- */
		bestwindow = GetBestWindow( from1, to1, bestdiagonal, sequence, lsequence, reducedots );
		if ( bestwindow.from == 0 && verbose > LL1) 
			printf ("-----> Iterating, because no windows found on best diagonal\n") ;

		lrepeat = bestwindow.to - bestwindow.from + 1;

		if ( lrepeat < minlrepeat && verbose > LL1) 
			printf("-----> Iterating, because window is too small (%i/%i)\n", lrepeat, minlrepeat);

		level_diagonal++;
	}

	if (bestwindow.from == 0) {
		if (verbose > LL1) 
			printf("-----> Leaving, because no best window found\n");
		return NULL;
	}

	/* align with wrapping -------------------------------------- */

	aali     = AlignBestWindow( bestwindow, lsequence, origdots );
	if ( aali->score < minswscore ) { 
		if (verbose > LL1) 
			printf("---> Leaving, because score is too low (%5.2f/%5.2f)\n", aali->score, minswscore);
		FreeArrayAliMemory( aali );
		return NULL;
	} 

	if (verbose > LL2) 
		DumpArrayAlignment (aali,"After Dotalign");

	/* convert to profile alignment ------------------------------------------- */
	MapAlignment( bestwindow.from, bestwindow.from, bestwindow.to, aali, 1);
	if (verbose > LL3) 
		DumpArrayAlignment( aali, "after remapping");

	/* eliminate unaligned profile positions, which occur, if parts of bestwindow are masked !!! */
	lrepeat = ReduceProfileAlignment ( aali, lrepeat, 1 );
	if (verbose > LL2) 
		printf("repeat length after reduction: %i\n", lrepeat);

	if ( lrepeat < minlrepeat) {
		if (verbose > LL2) printf("---> Leaving after reduction, because window is too small (%i/%i)\n", lrepeat, minlrepeat);
		for (i = 1; i < lsequence; i++) aali->ali[i] = 0;
		j = 1; for (i = bestwindow.from; i <= bestwindow.to; i++) aali->ali[i] = j++;
		*returnnrepeats = 1; 
		*returnscore    = 0;
		return aali;
	}

	/* Score alignment against itself */
	AutoScore( aali, lrepeat, 0, sequence ); 
	if (verbose > LL2) 
		DumpArrayAlignment( aali, "after reduction and rescoring");

	/* calculate the register  ----------------------------------------------------------*/
	NewestRegister ( aali, lrepeat, bestwindow.from ); 
	if (verbose > LL2) 
		DumpArrayAlignment (aali,"After Register");

	/* collect additional repeats -------------------------------------------*/
	niterations = 0;
	tempdots = CopyDotsHalf( origdots, UPPER );

	while (niterations < MAX_NITERATIONS) {
		niterations++;

		if (verbose > LL2) 
			printf("--->Iteration: %i\n", niterations);

		/* remap and mask dots 
	   in row: sequence
	   in col: profile
		 */
		tempdots2 = FoldDots( tempdots, origcounts, sequence, lsequence, aali, lrepeat, origprofile ) ; 

		if (verbose > LL3) {
			sprintf(buffer, "folddots_%i_%i_%i.png", iteration, level + 1, niterations );
		}

		if (verbose > LL4) {
			printf("Restrained dots\n");
			PrintDots( tempdots2->row, tempdots2->col, tempdots2->score, tempdots2->ndots);
		}

		if (tempdots2->ndots < minndots ) {
			if (verbose > LL1) 
				printf("---> Stopping repeat collection, because not enough dots (%i/%i)\n", tempdots2->ndots, minndots);
			FreeDotsMemory( tempdots2 );
			break;
		}

		/* align profile to sequence with wrapping */
		ali = AllocateAliMemory( MAXALI_SIZE);
		dotalign(   tempdots2->row, tempdots2->col, tempdots2->score, tempdots2->ndotsdiagonal, 
				1,			                   /* disallow rolling */
				1, lsequence,                                          
				1, lrepeat, 
				hp_gapopen_i, hp_gapelon_i,
				hp_gapopen_j, hp_gapelon_j, 
				ali );
		FreeDotsMemory( tempdots2 );


		if (verbose > LL3) {
			printf("---> dotalign: %5.2f\n", ali->score);
			DumpAlignment( ali, "Collecting-Alignment");
		}

		if (ali->score < mincollectscore) {
			if (verbose > LL1) 
				printf("---> Stopping repeat collection, because score is too low (%5.2f/%5.2f)\n", ali->score, mincollectscore);
			FreeAliMemory( ali );
			break;

		} else {
			if (verbose > LL1) 
				printf("---> Adding repeat (%i-%i) with length (%i) and score %5.2f\n",
						ali->align_j[ali->lastindex], ali->align_j[0], 
						ali->length, ali->score);

		}

		/* add newalignment to first alignment */
		for (i = 0; i <= ali->lastindex; i++) {
			aali->ali[ali->align_j[i]]    = ali->align_i[i];
			aali->scores[ali->align_j[i]] = ali->align_s[i];
		}

		if (verbose > LL2) 
			DumpArrayAlignment( aali, "added alignment");


		FreeAliMemory( ali );
	}

	FreeDotsMemory( tempdots );


	/* Score alignment against itself */
	AutoScore( aali, lrepeat, 0, sequence ); 

	if (verbose > LL1) { 
		printf("before register\n");
		PrintPrettyAli2( aali, sequence, lsequence);
	} 
	/*    PushGapVote ( bali, bestwindow.from, bestwindow.to, lsequence, repeatscores, &nrepeats ); */
	NewestRegister ( aali, lrepeat, bestwindow.from ); 

	if (verbose > LL1) {
		printf("after register\n");
		PrintPrettyAli2( aali, sequence, lsequence);
	}

	if (verbose > LL2) 
		DumpArrayAlignment( aali, "After repeat collection");

	/* FreeArrayAliMemory( bali ); */

	if (verbose > LL2) 
		DumpArrayAlignment( aali, "Alignment before pruning" );


	/* clean up alignment -----------------------------------------------------------*/
	repeatscores = Malloc( MAXNREPEATS * sizeof(CALCTYPE));
	zscores      = Malloc( MAXNREPEATS * sizeof(CALCTYPE));
	/*    intermediate_ali = PruneMultipleAlignment (aali, sequence, lrepeat, rg_gapopen, rg_gapelon, repeatscores, &nrepeats , &lrepeat_new);  */
	/*     intermediate_ali = SelectRepeats (aali,  */
	/* 				      sequence, lrepeat, */
	/* 				      rg_gapopen, rg_gapelon, */
	/* 				      repeatscores, &nrepeats,  */
	/* 				      &lrepeat_new ); */

	intermediate_ali = SelectRepeatsStatistical (
			aali, 
			sequence, lsequence, lrepeat,
			rg_gapopen, rg_gapelon,
			repeatscores, 
			zscores,
			&nrepeats, 
			&lrepeat_new 
	);
	FreeArrayAliMemory(aali);

	if (nrepeats <= 1) {
		if (verbose > 1) 
			printf("Leaving, because no/one repeat(s) left after pruning\n");
		nrepeats = 1;
		sprintf(resultbuffer, "");
	} else {

		if (verbose > LL2) {
			DumpArrayAlignment( intermediate_ali, "Final alignment" );
			printf ("Scores of %i repeats:\n", nrepeats);
			for (i = 0; i < nrepeats; i++) printf("%5.2f ", repeatscores[i]);
			printf("\n");
		}

		/* write alignment -------------------------------*/ 

		sumscore = 0; for (i=0; i < nrepeats; i++) sumscore += repeatscores[i];
		sprintf(resultbuffer, "---------------------------------------------------------------------------\n");
		sprintf(&resultbuffer[strlen(resultbuffer)], 
				"No. of Repeats|Total Score|Length  |Diagonal| BW-From|   BW-To|   Level\n");
		sprintf(&resultbuffer[strlen(resultbuffer)], 
				"%14i|%11.2f|%8i|%8i|%8i|%8i|%8i\n", 
				nrepeats, sumscore, lrepeat_new,
				bestdiagonal, bestwindow.from, bestwindow.to, iteration);
		sprintf(&resultbuffer[strlen(resultbuffer)], 
				"---------------------------------------------------------------------------\n");
		PrintPrettyMA( intermediate_ali, nrepeats, repeatscores, zscores, sequence, lsequence, 
				&resultbuffer[strlen(resultbuffer)] );
		sprintf(&resultbuffer[strlen(resultbuffer)], 
				"---------------------------------------------------------------------------\n");

		if (verbose > LL0) 
			printf(resultbuffer );

		/* check validity of the alignment -------------------------------------------------------------------- */

		/* if there are negative repeat-scores (fragmented alignment), use bestwindow for reduction */
		for (i = 0; i < nrepeats; i++) {
			if (repeatscores[i] < 0) 
				nrepeats = 1;
		}

	}

	free( repeatscores);

	/* do reduction -------------------------------------------------------------------------------- */
	/* if there is just one repeat, use bestwindow for reduction. Sometimes even the bestwindow does not align to
       itself, but is cut into several small pieces (disactivate this and look at pgbm_human */ 
	if (nrepeats == 1 ) {
		if (verbose > LL1) 
			printf("Using complete bestwindow for diagonal reduction\n");
		for (i = 1; i < lsequence; i++) intermediate_ali->ali[i] = 0;
		j = 1; for (i = bestwindow.from; i <= bestwindow.to; i++) intermediate_ali->ali[i] = j++;
		lrepeat_new   = bestwindow.to - bestwindow.from + 1;
		/* just copy the dots when looking for diagonals in smaller diagonals, excluding dots outside the window, because otherwise we have an infinte loop */
		tempdots = CopyDotsWindow( origdots, bestwindow.from + 1, bestwindow.to - 1 ); 
		/* discard = 1*/
	} else {
		/* remap dots, if we have a MA with several rows */
		/* tempdots = CreateEnhancedDots2( origdots, intermediate_ali, sequence, lsequence, lrepeat_new, origcounts, origprofile ) ; no remapping for titin-test*/
		tempdots = CopyDotsWindow( origdots, bestwindow.from + 1, bestwindow.to - 1 ); 
		/* look at best window */
	}
	free( zscores );

	smaller_alignment = Malloc( sizeof(char) * MAXRESULTSIZE );
	/* get aligned region */
	from2 = 0;         while ((from2 <= lsequence) && (intermediate_ali->ali[from2] == 0)) from2++;
	to2   = lsequence; while ((to2   > 0         ) && (intermediate_ali->ali[to2]   == 0)) to2--;

	/* call self again */
	smaller_ali = CalculateRepeatMA( from2, to2,
			sequence, lsequence,
			origprofile, lprofile,
			origdots, origcounts,
			level+1,
			&smaller_score,
			&smaller_nrepeats,
			smaller_alignment,
			tempdots);
	FreeDotsMemory ( tempdots );

	/* decide, which alignment to return */
	if (smaller_ali != NULL &&						/* return smaller alignment, if there is an alignment */
			(smaller_nrepeats >= nrepeats) 					/* and I found more repeats */
			/*(CompareAlignments( smaller_ali, intermediate_ali) ||                          */
			/* or this alignment would be discarded anyway */
	) {  
		*returnscore    = smaller_score; 
		*returnnrepeats = smaller_nrepeats;
		strcpy( returnalignment, smaller_alignment);
		free( smaller_alignment );
		FreeArrayAliMemory( intermediate_ali);
		return smaller_ali;
	} else {								/* return this alignment, which means defaulting to bestwindow */
		*returnscore    = sumscore;
		*returnnrepeats = nrepeats;
		strcpy( returnalignment, resultbuffer );
		free( smaller_alignment);
		if (smaller_ali != NULL) 
			FreeArrayAliMemory( smaller_ali);
		return intermediate_ali;
	}						       

}

/* run radar 
 * 
 * warning: dots are modified here - a diagonal is added.
 * */
int run_radar( 
			int lsequence,
			SEQTYPE * sequence,
			int lprofile,
			COUNTCOLUMN * counts,
			PROFILECOLUMN * profile,
			DOTS * dots )
{
	SEQTYPE * masked_sequence;
	DOTS *workdots;
	int from1, to1;
	int nrepeats, i;
	CALCTYPE sumscore;
	int masked;
	int nextlevel, level;
	char buffer[100];
	ARRAYALI *final_ali;
	FREQUENCYCOLUMN *frequencies;
	
	// not enough dots - no result
	if (dots->ndots < minndots) 
	{
		printf("No repeats found\n");
		return 0;
	}

	/* allocate memory for masked sequence */
	masked_sequence = (SEQTYPE*)Malloc( ( lsequence + 1 ) * sizeof(char));
	
	/* add diagonal */
	frequencies = AllocateFrequenciesMemory( lprofile );
	Counts2Frequencies( counts, lprofile, frequencies );
	AddDiagonal( dots, profile, frequencies, lprofile, sequence, lsequence );
	FreeFrequenciesMemory( frequencies );

	iteration     = 0;  /* starting with diagonal iteration + 1 */
	level         = 0;   
	nextlevel     = 1;  /* enter into first iteration right away */

	/* ---- start of main loop ---------------------------------------------------*/
	if (verbose > LL1) 
		printf( "--> Start of main loop\n" );
	
	while(1) 
	{

		iteration++;
		if (verbose > LL1) 
			printf ("--->Iterating: %i <-------------------------------------------------------------\n", iteration);

		/* enter a new iteration ---------------------------------------------- */
		if (nextlevel) 
		{
			level ++;
			if ( level > maxlevels) break;
			if (verbose > LL1) 
				printf ("--->Iterating: %i <-------------------------------------------------------------\n", level);
			memcpy( masked_sequence, sequence, sizeof(SEQTYPE) *(lsequence + 1)); /* unmask sequence */
			nextlevel = 0;
		}

		workdots = MaskDots(dots, masked_sequence, lsequence );
		if (workdots->ndots < minndots) 
		{
			if (verbose > LL2) 
				printf ("Exiting, because not enough dots\n");
			nextlevel = 1;
			continue;
		}		

		if (verbose > LL4)
			PrintDots( workdots->row, workdots->col, workdots->score, workdots->ndotsdiagonal);

		if (verbose > LL3) 
		{
			sprintf(buffer, "initialdots_%i.png", iteration);
		}

		if (verbose > LL1) 
			printf("Using %i (%i) dots\n", workdots->ndots, workdots->ndotsdiagonal);

		from1 = 1;
		to1   = lsequence;

		/* start repeatfiextern int iteration;
extern char globalresult[MAXRESULTSIZE];
		 * nding ------------------------------------------------- */
		final_ali = CalculateRepeatMA( 
				from1, to1, 
				sequence, lsequence,
				profile, lprofile,
				workdots, counts,
				0,
				&sumscore,
				&nrepeats,
				&globalresult[strlen(globalresult)],
				workdots);
		
		FreeDotsMemory( workdots );

		if (final_ali == NULL) 
		{
			nextlevel = 1;
			continue;
		};

		if (verbose > LL0) 
			printf("Alignment with score %5.2f and %i repeats used for masking\n", sumscore, nrepeats);

		/* mask dots and iterate */
		masked = 0;
		for (i = 1; i <= final_ali->length; i++) 
		{
			if (final_ali->ali[i] != 0) 
			{
				masked_sequence[i] = MASKCODE;
				masked = 1;
			}
		}
		FreeArrayAliMemory( final_ali );

		if (!masked) 
		{
			if (verbose > LL2) 
				printf("Nothing masked, exiting level\n");
			nextlevel = 1;
		}
	}

	free( masked_sequence);

	return 1;
}

int radar_run_from_files( 
			 const char * filename_sequence,
			 const char * filename_ma,
			 const char * filename_lfasta,
			 const char * filename_lfasta2,
			 unsigned int random_seed)
{
	int lprofile;
	PROFILECOLUMN *profile;
	COUNTCOLUMN   *counts;
	SEQTYPE * sequence;
	int lsequence;
	DOTS *dots;
	char buffer[100];

	if (verbose > LL1)
		printf("Starting RADAR Version $Name:  $\n");		

	WellcomeMsg();
	globalresult[0] = '\0'; 

	if (random_seed != 0) 
	    srand(random_seed);
	else
	    srand((unsigned int)clock());

	/* read sequence and allocate memory for the masked sequence -----*/
	if (verbose > LL1) printf("Reading sequence...");
	lsequence = ReadSequence( filename_sequence, &sequence );
	if (verbose > LL1) printf("Done (length=%i)\n", lsequence);
	if (verbose > LL3) 
	{
		printf ("read sequence of length %i:\n", lsequence);
		PrintSequence( sequence, lsequence );
	}

	counts  = AllocateCountsMemory( lsequence );
	profile = AllocateProfileMemory( lsequence );

	if (verbose > LL1) printf("Reading counts from multiple alignment...");
	ReadCountsFromMA( filename_ma, counts ); 
	if (verbose > LL1) printf("Done\n"); 

	if (verbose > LL1) printf("Calculating profile...");
	Counts2Profile( counts, profile, lsequence);
	if (verbose > LL1) printf("Done\n");

	lprofile = lsequence;

	if (verbose > LL3) 
	{
		PrintCounts ( counts,  lprofile );
		PrintProfile( profile, lprofile);
	}

	/* free( ma ); */
	/* read lfasta-file and create dotsfile in dots ------------------------------------------------------------------- */
	dots = AllocateDotsMemory( MAXNDOTS );
	if (verbose > LL1) printf("Reading dots...");
	if (strlen(filename_lfasta2) > 0)
		dots->ndots = FastaFiles2DotsCorrectSaveSpace( dots->row,
				dots->col, 
				dots->score,  
				filename_lfasta, 
				filename_lfasta2, 
				profile, 
				counts, 
				lprofile, 
				mindiag);
	else
		dots->ndots = Fasta2DotsCorrect( dots->row, 
				dots->col, 
				dots->score, 
				filename_lfasta,
				profile, 
				counts, 
				lprofile, 
				mindiag);

	if (verbose > LL1) 
		printf("Done\n");

	if (verbose > LL3)
		PrintDots( dots->row, dots->col, dots->score, dots->ndotsdiagonal);
	

	// run radar
	int success = run_radar( 
			lsequence, 
			sequence,
			lprofile,
			counts,
			profile,
			dots );
	
	// output final result
	if (verbose > LL1) 
		printf("---Final result----------------------------------------------------------\n");

	if (strlen( globalresult ) > 0)
		printf(globalresult);
	else
		printf("No repeats found\n");
	/* printf("radar finished without problems\n"); web-version*/

	free( sequence );
	FreeCountsMemory( counts );
	FreeProfileMemory( profile );
	FreeDotsMemory( dots );
	
	return EXIT_SUCCESS;
}

