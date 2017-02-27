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
$Id: repeatstatistics.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#include "radar.h"

void CalculateMeanAndStdDev( scores, size, mean, standard_deviation) 
CALCTYPE *scores;
int       size;
double   *mean;
double   *standard_deviation;
{
    int i;
    double total;
    double m, s;
    
    /* calculate mean */
    total = 0; for (i = 0; i < size; i++) total+= (double)scores[i];
    m = total / (double)size;

    /* calculate standard-deviation */
    total = 0; for (i = 0; i < size; i++) total+= ((double)scores[i] - m) * ((double)scores[i] - m);
    s = sqrt( total / (double)size);
    *mean = m;
    *standard_deviation = s;
}

int Rand ( max ) 
    int max;
{
    return (1+(int)((double)max*rand()/(RAND_MAX+1.0)));
}

 

void ShuffleSequence( rsequence, lsequence ) 
    SEQTYPE *rsequence;
int lsequence;
{
    int i,j; 
    SEQTYPE x;
    
    for (i = lsequence; i >= 1; i--) {
	j = Rand(lsequence);
	x = rsequence[j];
	rsequence[j] = rsequence[i];
	rsequence[i] = x;
    }
}    


int * CalculateCounts( sequence, lsequence )
SEQTYPE * sequence;
int lsequence;
{
    int i;
    int *counts = Malloc ( sizeof (int ) * PROFILEWIDTH );
    for (i = 0; i < PROFILEWIDTH; i++) counts[i] = 0;
    for (i = 1; i <= lsequence; i++)   counts[sequence[i]] ++;

    return counts;
}

double * __Counts2Frequencies( counts )
int * counts;
{
    int i, total = 0;
    double *frequencies = Malloc ( sizeof (double) * PROFILEWIDTH );
    
    for (i = 0; i < PROFILEWIDTH; i++) total += counts[i];
    for (i = 0; i < PROFILEWIDTH; i++) frequencies[i] = (double)counts[i] / (double)total;
    
    return frequencies;
}


void CalculateZScore(  
		     sequence, lsequence, 
		     profile, lprofile, 
		     gapopen_i, gapelon_i,
		     gapopen_j, gapelon_j,
		     iterations,
		     mean, standard_deviation
		     ) 
    SEQTYPE * sequence;
int lsequence;
PROFILECOLUMN *profile;
int lprofile;
CALCTYPE gapopen_i, gapelon_i;
CALCTYPE gapopen_j, gapelon_j;
int iterations;
double *mean, *standard_deviation;
{
    SEQTYPE  *rsequence   = CopySequence( sequence, lsequence );
    CALCTYPE *scores      = Malloc( iterations * sizeof(CALCTYPE));
    int i;

    /* get scores for iterations alignments of shuffled sequences against profile---------------------------------------------------*/
    for (i = 0; i < iterations; i++) {
	ShuffleSequence( rsequence, lsequence );
	scores[i] = lps_align_score( profile, rsequence, 
				     lprofile, lsequence, 
				     -gapopen_i, -gapelon_i, 
				     -gapopen_j, -gapelon_j);
	/* printf("Alignment score %5.2f\n", scores[i]); */
    }
    
    /* calculate mean and standard-deviation---------------------------------------------------------------------------------------- */
    CalculateMeanAndStdDev( scores, iterations, mean, standard_deviation);
    free (scores);
    FreeSequenceMemory (rsequence);
    
}
