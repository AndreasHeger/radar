/*
RADAR - a program to rapidly detect and align repeats in protein sequences
Copyright (C) 2000 Andreas Heger

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

##########################################################################
$Id: getmaxdiag.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/
#include "align.h"
#include "align_tools.h"
#include "toolbox.h"
#include <math.h>


int getnbestdiagpostrace( row, col, score, ndots, mincol, maxcol, minrow, maxrow, level, region, minscore, bestscore)
int     *col;
int     *row;
CALCTYPE *score;
int      ndots;
int      mincol, maxcol, minrow, maxrow; 
int      level;
int region;
CALCTYPE minscore;
CALCTYPE *bestscore;

{ 

    int i;
    CALCTYPE * count;
    CALCTYPE * bestc;
    int      * index;
    int d, d1, n;
    int maxdiag, maxindex;
    int bestdiag;
    int dmax;
    
    maxdiag = ((maxrow - mincol) > (maxcol - mincol)) ? (maxrow - mincol) : (maxcol - minrow);
    
    if (maxdiag < 0 ) maxdiag = -maxdiag;
#ifdef DEBUG
    printf("GetNBestDiagTrace called with mincol %i, maxcol %i , minrow %i, maxrow %i, ndots %i level %i, region %i, minscore %5.2f\n",
	   mincol, maxcol, minrow, maxrow, ndots, level, region, minscore) ;
    printf("Maximum Diagonal: %i\n", maxdiag);
    PrintDots(row, col, score, ndots); 
#endif

    maxindex = maxdiag;
    
    count = (CALCTYPE*)Malloc( (maxindex + 1) * sizeof(CALCTYPE) ); /* allocate memory for current scores */
    bestc = (CALCTYPE*)Malloc( (maxindex + 1) * sizeof(CALCTYPE) ); /* allocate memory for best scores */
    index = (int*)     Malloc( (maxindex + 1) * sizeof(int) );      /* allocate memory for index to sort */
    
    for (i = 0; i <= maxindex; i++) count[i] = 0;  
    for (i = 0; i <= maxindex; i++) bestc[i] = 0;  
    for (i = 0; i <= maxindex; i++) index[i] = i;  

    for (i = 0; i < ndots; i++) {
	if (row[i] > col[i]) { 
	    continue;
	}
	if (row[i] >= minrow &&
	    row[i] <= maxrow &&
	    col[i] >= mincol &&
	    col[i] <= maxcol ) {
	    d = col[i] - row[i];
	    if (score[i] > 0 ) {
		count[d] += score[i];            /* add positive scores to current scores */
	    }	    else 
		if (count[d] > bestc[d]) {
		    bestc[d] = count[d];          /* end summing up, if negative residue starts */
		    count[d] =0;                  /* reset trace to start from zero */
		}
	}
    }

    /* transfer last values */
    for (i = 0; i < maxindex; i++) 
	if (bestc[i] < count[i])
	    bestc[i] = count[i];
    
    /* add up over region, note due to this calculation the index of the center shifts by one back, i.e. in 
     in row 31 is the added up value for row 32 */

    region += region;
    for (d = 0; d <= maxindex; d++) {
	dmax = ((d + region) >= maxindex ) ? maxindex : (d + region);
	for (d1 = d + 1; d1 <= dmax; d1++ )    /* sum over many diagonals */
	    bestc[d] += bestc[d1];
    }
    
    /* sort diagonals */
    QuickSortIndicesDouble( 0, maxindex - 1, index, bestc);
#ifdef DEBUG  
    printf("Sorted index\n");
    for (i = 0; i < maxindex; i++) {
	printf("%i %i %5.2f\n", i, index[i], bestc[index[i]]);
    }
#endif 
    
    /* if there are several diagonals with the same score, pick the middle one */
    bestdiag   = index[maxindex - level - 1 ];
    /* calculate average diagonal of diagonals with same score around best diagonal */
    d = bestdiag;
    n = 1;
    i = bestdiag + 1;
    while ( (i < maxindex) && (bestc[i] == bestc[bestdiag])){ d += i; n++; i++;}
    i = bestdiag - 1; 
    while ( (i >= 0 ) &&      (bestc[i] == bestc[bestdiag])){ d += i; n++; i--;}
    
    bestdiag   = rint(d / n); 
    *bestscore = bestc[bestdiag];
    bestdiag  -= 1; /* correction for shift above */	
    free ( bestc );
    free ( count );
    free ( index );
    return bestdiag;
}




/*--------------------------------------------------------------------
  Imagine the following dotplot:

  12345678
 1\       A
 2_\______
 3\ \     B
 4_\_\____
 5\ \ \   C
 6_\_\_\___
 7\ \ \ \ D
 8_\_\_\_\_
    a b c d
 with the four repeats A, B, C, D with a length of 2 each.
 if you call getnbestdiag with from = 1 and to = 8 and level =
 2, you retrieve the offset of diagonal b.
 If you want to pick up diagonal a, either look for the third highest
 diagonal from d or for the highest scoring diagonal from b, i.e. looking
 in the box (1,6);(3,8). 
 The algorithm just looks at diagonals below the main diagonal, i.e. row has 
 is smaller as col. Remember this when specifying mincol, minrow, etc.
*/

int getnbestdiagsum( row, col, score, ndots, mincol, maxcol, minrow, maxrow, level, region, minscore, bestscore)
int     *col;
int     *row;
CALCTYPE *score;
int      ndots;
int      mincol, maxcol, minrow, maxrow; 
int      level;
int region;
CALCTYPE minscore;
CALCTYPE *bestscore;

{ 

    int i;
    CALCTYPE * count;
    int      * index;
    int d, d1, n;
    int maxdiag, maxindex;
    int bestdiag;
    int dmin, dmax;
    
    maxdiag = ((maxrow - mincol) > (maxcol - mincol)) ? (maxrow - mincol) : (maxcol - minrow);
    
    if (maxdiag < 0 ) maxdiag = -maxdiag;
#ifdef DEBUG
    printf("GetNBestDiagSum called with mincol %i, maxcol %i , minrow %i, maxrow %i, ndots %i level %i, region %i, minscore %5.2f\n",
	   mincol, maxcol, minrow, maxrow, ndots, level, region, minscore) ;
    printf("Maximum Diagonal: %i\n", maxdiag);
    PrintDots(row, col, score, ndots); 
#endif

    maxindex = maxdiag;

    count = (CALCTYPE*)Malloc( (maxindex + 1)* sizeof(CALCTYPE) ); /* allocate memory */
    index = (int*)     Malloc( (maxindex + 1)* sizeof(int) );      /* allocate memory */

    for (i = 0; i < maxindex; i++) count[i] = 0;  
    for (i = 0; i < maxindex; i++) index[i] = i;  

    for (i = 0; i < ndots; i++) {
	if (row[i] > col[i]) { 
	    continue;
	}
	if (row[i] >= minrow &&
	    row[i] <= maxrow &&
	    col[i] >= mincol &&
	    col[i] <= maxcol ) {
	    d = col[i] - row[i];
	    dmin = ((d - region) < 0)          ?        0 : (d - region);
	    dmax = ((d + region) >= maxindex ) ? maxindex : (d + region);
	    for (d1 = dmin; d1 <= dmax; d1++ ) /* sum over many diagonals */
		if (score[i] > minscore)       /* only use scores above threshold, because bad traces otherwise mask out good diagonal */
		    count[d1] += score[i];
	}
    }
    
    QuickSortIndicesDouble( 0, maxindex - 1, index, count);
#ifdef DEBUG  
    printf("Sorted index\n");
    for (i = 0; i < maxindex; i++) {
	printf("%i %i %5.2f\n", i, index[i], count[index[i]]);
    }
#endif 

    /* if there are several diagonals with the same score, pick the middle one */
    bestdiag   = index[maxindex - level - 1 ];
    /* calculate average diagonal of diagonals with same score around best diagonal */
    d = bestdiag;
    n = 1;
    i = bestdiag + 1;
    while ( (i < maxindex) && (count[i] == count[bestdiag])){ d += i; n++; i++;}
    i = bestdiag - 1;
    while ( (i >= 0 ) &&      (count[i] == count[bestdiag])){ d += i; n++; i--;}
    
    bestdiag   = rint(d / n);	
    *bestscore = count[bestdiag];
    free ( count );
    free ( index );
    return bestdiag;
}


int getnbestdiag( row, col, score, ndots, mincol, maxcol, minrow, maxrow, level, bestscore)
int     *col;
int     *row;
CALCTYPE *score;
int      ndots;
int      mincol, maxcol, minrow, maxrow; 
int      level;
CALCTYPE *bestscore;
{ 

    int i;
    CALCTYPE * count;
    int      * index;
    int d;
    int maxdiag, maxindex;
    int bestdiag;
    
    maxdiag = ((maxrow - mincol) > (maxcol - mincol)) ? (maxrow - mincol) : (maxcol - minrow);
    
    if (maxdiag < 0 ) maxdiag = -maxdiag;

#ifdef DEBUG
    printf("getnbestdiag called with mincol %i, maxcol %i , minrow %i, maxrow %i, ndots %i, level %i\n",
	   mincol, maxcol, minrow, maxrow, ndots, level) ;
    printf("Maximum Diagonal: %i\n", maxdiag);
    PrintDots(row, col, score, ndots); 
#endif

    maxindex = maxdiag;

    count = Malloc( (maxindex + 1) * sizeof(CALCTYPE) ); /* allocate memory */
    index = Malloc( (maxindex + 1) * sizeof(int) );      /* allocate memory */
    
    for (i = 0; i < maxindex; i++) count[i] = 0;  
    for (i = 0; i < maxindex; i++) index[i] = i;  
    
    for (i = 0; i < ndots; i++) {
	if (row[i] > col[i]) { 
	    continue;
	}
	if (row[i] >= minrow &&
	    row[i] <= maxrow &&
	    col[i] >= mincol &&
	    col[i] <= maxcol ) {
	    d = col[i] - row[i];
	    count[d] += score[i];
	}
    }

    QuickSortIndicesDouble( 0, maxindex - 1, index, count);

#ifdef DEBUG
    printf("Sorted index\n");
    for (i = 0; i < maxindex; i++) {
	printf("%i %i %5.2f\n", i, index[i], count[index[i]]);
    }
#endif
    bestdiag   = index[maxindex - level - 1];
    *bestscore = count[bestdiag];
    free ( count );
    free ( index );
    return bestdiag;
}


int getmaxdiag ( row, col, score, ndots, mincol, maxcol, minrow, maxrow, bestscore)
int     *col;
int     *row;
CALCTYPE *score;
int      ndots;
int      mincol, maxcol, minrow, maxrow; 
CALCTYPE *bestscore;
{ 

    int i;
    CALCTYPE * count;
    int d;
    int maxdiag, maxindex;
    int best, bestdiag;
    
    maxdiag = ((maxrow - mincol) > (maxcol - mincol)) ? (maxrow - mincol) : (maxcol - minrow);
    
    if (maxdiag < 0 ) maxdiag = -maxdiag;

#ifdef DEBUG
    printf("GetMaxDiag called with mincol %i, maxcol %i , minrow %i, maxrow %i, ndots %i\n",
	   mincol, maxcol, minrow, maxrow, ndots) ;
    printf("Maximum Diagonal: %i\n", maxdiag);
    PrintDots(row, col, score, ndots); 
#endif

    maxindex = 2 * maxdiag + 1;

    count = Malloc( maxindex * sizeof(CALCTYPE) ); /* allocate memory */

    for (i = 0; i < maxindex; i++) count[i] = 0;  

    for (i = 0; i < ndots; i++) {
	if (row[i] > col[i]) { 
	    continue;
	}
	if (row[i] >= minrow &&
	    row[i] <= maxrow &&
	    col[i] >= mincol &&
	    col[i] <= maxcol ) {
	    d = row[i] - col[i] + maxdiag;
	    count[d] += score[i];
	}
    }
    
    best     = 0;
    bestdiag = 0;
    for (i = 0; i < maxindex; i++) {
	if (count[i] > best) {
	    best     = count[i];
	    bestdiag = i;
	    *bestscore = count[i];
	}
    }
    
    bestdiag -= maxdiag;
	
    free ( count );
    return bestdiag;
}
