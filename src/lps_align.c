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
$Id: lps_align.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#include "align.h"
#include "align_tools.h"
#include "toolbox.h"

static CALCTYPE gi, hi, mi;				/* g = G, h = H, m = g+h */
static CALCTYPE gj, hj, mj;				/* g = G, h = H, m = g+h */

static CALCTYPE CC[NMAX+1], DD[NMAX+1];	/* Forward cost-only vectors */

static CELL H1[NMAX+1];			/* Pointer Array for start of local alignment */
static CELL H2[NMAX+1];			/* Pointer Array for start of local alignment */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */


static CALCTYPE __lps_align( A,B,M,N,tb,te, startcell, stopcell) 
PROFILECOLUMN *A; 
SEQTYPE       *B; 
int            M, N;
CALCTYPE      tb, te; 
CELL *startcell, *stopcell;

{ register int   i, j;
  CALCTYPE c, e, d, s;
  CALCTYPE *wa;
  CALCTYPE max;

	   max = 0;
	   CC[0] = 0;
	   H1[0].j = 1;
	   H1[0].i = 1;
	   for (j = 1; j <= N; j++) { 
	       CC[j] = 0;
	       DD[j] = -gi;
	       H2[j].i = 1;
	       H2[j].j = j + 1;
	   }

	   for (i = 1; i <= M; i++) {
	       H1[0].i = i;
	       H1[0].j = 1;
	       s = CC[0];
	       CC[0] = c = 0;
	       e = -gj;
	       wa = A[i];
	       for (j = 1; j <= N; j++) { 

		   if ((c =   c   - mj) > (e =   e   - hj)) e = c;
		   if ((c = CC[j] - mi) > (d = DD[j] - hi)) d = c;
		   c = s + wa[B[j]];
		   H1[j].i = H2[j].i;
		   H1[j].j = H2[j].j;
		   if (e > c) c = e;
		   if (d > c) c = d;
		   if (c <= 0) {
		       c = 0; 		/* the local alignment part */
		       H2[j].i = i + 1;	/* if new alignment starts, set pointer to current cell if matched */
		       H2[j].j = j + 1; 
#ifdef DEBUG
	      	       printf("|o"); 
#endif
		   } else {
		       if ( c == e ) {
			   H2[j].i = H2[j-1].i; /* propagate pointer for deletion*/ 
			   H2[j].j = H2[j-1].j; 
#ifdef DEBUG
			   printf("|i"); 
#endif
		       } else if ( c == d ) {
			   H2[j].i = H1[j].i;    /* propagate pointer for insertion*/ 
			   H2[j].j = H1[j].j; 
#ifdef DEBUG
			   printf("|d"); 
#endif
		       } else {
			   H2[j].i = H1[j-1].i;     /* propagate pointer for match*/ 
			   H2[j].j = H1[j-1].j; 
#ifdef DEBUG
			   printf("|m"); 
#endif
		       }
		   }
		   s = CC[j];
		   CC[j] = c;

#ifdef DEBUG
		   printf("%5.2f %5.2f", c, CC[j]); /*  , H2[j].i, H2[j].j); */
#endif
		   if (max < c) { /* save max */  
		       max = c;
 		       stopcell->i = i;  
 		       stopcell->j = j;  
		       startcell->i = H2[j].i; 
		       startcell->j = H2[j].j; 
		   }
		   DD[j] = d;
	       }
#ifdef DEBUG
	       printf ("\n");
#endif
	   }
	   
	   /* startcell and stopcell contain the boundaries  cell of the local alignment */
	   return max;
}


void 
lps_align(A, B, M, N, gopi, gepi, gopj, gepj, ali)
PROFILECOLUMN *A;
SEQTYPE *B;
int      M, N;
CALCTYPE gopi, gepi, gopj, gepj;
ALI     *ali;
{

    int i;
    CELL startcell, stopcell;
    CALCTYPE score;
    /* Setup global parameters */

#ifdef DEBUG
    printf("lps_align called with gopi %5.2f, gepi %5.2f, gopj %5.2f, gepj %5.2f, M %i, N %i\n", gopi, gepi, gopj, gepj, M, N);
    PrintProfile( A, M);
    PrintSequence( B, N);
#endif

    gi = gopi; hi = gepi; mi = gi+hi;
    gj = gopj; hj = gepj; mj = gj+hj;

    /* do the local alignment to find start and end */
    score = __lps_align(A, B, M, N,
			-gi, -gj,
			&startcell, &stopcell);

#ifdef DEBUG
    printf("Startcell: i %i j %i; Stopcell: i %i j %i; Score %5.2f\n", 
	   startcell.i, startcell.j,
	   stopcell.i , stopcell.j, 
	   score);
#endif
    /* to get rid of all the -1s */
    startcell.i --;
    startcell.j --;
    
    /* do global alignment to find the trace between start and end */
    gps_align(
	      A + startcell.i,      /* sequence starts with a 1, array starts with 0 */
	      B + startcell.j,	/* thererfore subtract 1 so that there is one dummy*/
	      stopcell.i - startcell.i,
	      stopcell.j - startcell.j,
	      gopi, gepi, 
	      gopj, gepj, 
	      ali); 
    
    ali->score = score;

    for (i = 0; i <= ali->lastindex; i++) {
        ali->align_i[i] += startcell.i;
        ali->align_j[i] += startcell.j;
    }

}

CALCTYPE
lps_align_score(A, B, M, N, gopi, gepi, gopj, gepj, ali)
PROFILECOLUMN *A;
SEQTYPE *B;
int      M, N;
CALCTYPE gopi, gepi, gopj, gepj;
{
    CELL startcell, stopcell;
    CALCTYPE score;
    /* Setup global parameters */

#ifdef DEBUG
    printf("lps_align called with gopi %5.2f, gepi %5.2f, gopj %5.2f, gepj %5.2f, M %i, N %i\n", gopi, gepi, gopj, gepj, M, N);
    PrintProfile( A, M);
    PrintSequence( B, N);
#endif
    
    gi = gopi; hi = gepi; mi = gi+hi;
    gj = gopj; hj = gepj; mj = gj+hj;
    
    /* do the local alignment to find start and end */
    score = __lps_align(A, B, M, N,
			-gi, -gj,
			&startcell, &stopcell);
#ifdef DEBUG
    printf("Startcell: i %i j %i; Stopcell: i %i j %i; Score %5.2f\n", 
	   startcell.i, startcell.j,
	   stopcell.i , stopcell.j, 
	   score);
#endif

    return score;
}
