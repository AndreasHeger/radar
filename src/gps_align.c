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
$Id: gps_align.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

/** Part of this subroutine has been copied and modified from the FASTA package:

    W. R. Pearson and D. J. Lipman (1988),
    "Improved Tools for Biological Sequence Analysis", PNAS 85:2444-2448
    
    version 2.0u66 September 1998
    
    X. Huang and W. Miller (1991) Adv. Appl. Math. 12:373-381

    This is the excerpt from the COPYRIGHT-file:

    Copyright 1988, 1991, 1992, 1993, 1994 1995, by William
    R. Pearson and the University of Virginia.  All rights
    reserved. The FASTA program and documentation may not be sold or
    incorporated into a commercial product, in whole or in part,
    without written consent of William R. Pearson and the University
    of Virginia.  For further information regarding permission for
    use or reproduction, please contact:
 
    David Hudson
    Assistant Provost for Research
    University of Virginia
    P.O. Box 9025
    Charlottesville, VA  22906-9025
 
    (804) 924-6853
 
    COPYRIGHT line 1/17 (END)                   

*/

#include "align.h"
#include "align_tools.h"
#include "toolbox.h"

static CALCTYPE gi, hi, mi;				/* g = G, h = H, m = g+h */
static CALCTYPE gj, hj, mj;				/* g = G, h = H, m = g+h */

static int *sapp;				/* Current script append ptr */
static int  last;				/* Last script op appended */
						/* Append "Delete k" op */
static CALCTYPE CC[NMAX+1], DD[NMAX+1];	/* Forward cost-only vectors */
static CALCTYPE RR[NMAX+1], SS[NMAX+1];	/* Reverse cost-only vectors */

/* align(A,B,M,N,tb,te) returns the cost of an optimum conversion between
   A[1..M] and B[1..N] that begins(ends) with a delete if tb(te) is zero
   and appends such a conversion to the current script.                   */


/* 
   row: A, M, tb, i
   col: B, N, te, j
*/   

static CALCTYPE __gps_align(A,B,M,N,tb,te,topr,botr,lc,rc)  
PROFILECOLUMN *A; 
SEQTYPE       *B; 
int           M, N;
CALCTYPE      tb, te; 
char topr, botr, lc, rc;

{        
    int   midi, midj, type;	/* Midpoint, type, and cost */
    CALCTYPE midc;
    { 
	register int   i, j;
	CALCTYPE  c, e, d, s;
	CALCTYPE  t;
	CALCTYPE *wa;

	/* Boundary cases: M <= 1 or N == 0 */

	if (N <= 0) {                   /* if length of B is zero */
	    if (M > 0) DEL(M);          
	    if (topr || botr)           /* no gap penalties, if in top or bottom row */
		return 0;
	    else 
		return -gapi(M);
	}
	if (M <= 1) {                  
	    if (M <= 0) {              /* if length of A is zero, i.e. no rows */
		INS(N);
		if (topr || botr)
		    return 0;          /* no gap penalties, if in top or bottom row */
		else 
		    return -gapj(N);
	    }
	    if (topr) {                /* if length of A is one */
		if (rc)                /* if in toprow and rightcolumn */
		    midc = 0;
		else 
		    midc = te-hj;
		midj = 0;
		wa = A[1];
		for (j = 1; j <= N; j++) { 
		    c = wa[B[j]] - gapj(N-j);
		    if (c > midc) { 
			midc = c;
			midj = j;
		    }
		}
	    } else if (botr) {
		if (lc) 
		    midc = 0;
		else 
		    midc = tb-hj;
		midj = 0;
		wa = A[1];
		for (j = 1; j <= N; j++) { 
		    c = -gapj(j-1) + wa[B[j]];
		    if (c > midc) { 
			midc = c;
			midj = j;
		    }
		}
	    } else {
		if (tb < te) tb = te;
		if (lc || rc) 
		    midc = -gapj(N);
		else 
		    midc = (tb-hj) - gapj(N);
		midj = 0;
		wa = A[1];
		for (j = 1; j <= N; j++) { 
		    c = -gapj(j-1) + wa[B[j]] - gapj(N-j);
		    if (c > midc) {
			midc = c;
			midj = j;
		    }
		}
	    }
	    if (midj == 0) { 
		INS(N) DEL(1);
	    } else { 
		if (midj > 1) INS(midj-1);
		REP;
		if (midj < N) INS(N-midj);
	    }
	    return midc;
	}
	
/* Divide: Find optimum midpoint (midi,midj) of cost midc */

  midi = M/2;			/* Forward phase:                          */
  CC[0] = 0;			/*   Compute C(M/2,k) & D(M/2,k) for all k */
  if (topr) {
      for (j = 1; j <= N; j++) { 
	  CC[j] = 0;
	  DD[j] = -gi;
      }
  } else {
      t = -gj;
      for (j = 1; j <= N; j++) { 
	  CC[j] = t = t-hj;
       	  DD[j] = t-gi;
      }
  }
  t = tb;
  for (i = 1; i <= midi; i++) { 
      s = CC[0];
      if (lc) {
	  CC[0] = c = 0;
	  e = -gj;
      } else {
	  CC[0] = c = t = t-hj;
	  e = t-gj;
      }
      wa = A[i];
      for (j = 1; j <= N; j++)
        { if ((c =   c   - mj) > (e =   e   - hj)) e = c;
	  if ((j == N) && rc) {
             if ((c = CC[j]) > (d = DD[j])) d = c;
	  } else {   
             if ((c = CC[j] - mi) > (d = DD[j] - hi)) d = c;
	  }
          c = s + wa[B[j]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = CC[j];
          CC[j] = c;
          DD[j] = d;
        }
    }
  DD[0] = CC[0];

  RR[N] = 0;			/* Reverse phase:                          */
  				/*   Compute R(M/2,k) & S(M/2,k) for all k */
  if (botr) {
	for (j = N-1; j >= 0; j--)
	  { RR[j] = 0;
	    SS[j] = -gi;
          }
  } else {
      t = -gj;
      for (j = N-1; j >= 0; j--) { 
	  RR[j] = t = t-hj;
      	  SS[j] = t-gi;
      }
  }
  t = te;
  for (i = M-1; i >= midi; i--)
    { s = RR[N];
      if (rc) {
	RR[N] = c = 0;
	e = -gj;
      } else {
      	RR[N] = c = t = t-hj;
      	e = t-gj;
      }
      wa = A[i+1];
      for (j = N-1; j >= 0; j--)
        { if ((c =   c   - mj) > (e =   e   - hj)) e = c;
	  if ((j == 0) && lc) {
             if ((c = RR[j]) > (d = SS[j])) d = c;
	  } else {
             if ((c = RR[j] - mi) > (d = SS[j] - hi)) d = c;
	  }
          c = s + wa[B[j+1]];
          if (e > c) c = e;
          if (d > c) c = d;
          s = RR[j];
          RR[j] = c;
          SS[j] = d;
        }
    }
  SS[N] = RR[N];

  midc = CC[0]+RR[0];		/* Find optimal midpoint */
  midj = 0;
  type = 1;
  for (j = 0; j <= N; j++)
    if ((c = CC[j] + RR[j]) >= midc)
      if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
        { midc = c;
          midj = j;
        }
  if (rc) {
    if ((c = DD[N] + SS[N]) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  } else {
    if ((c = DD[N] + SS[N] + gi) > midc)
      { midc = c;
        midj = N;
        type = 2;
      }
  }
  for (j = N-1; j > 0; j--)
    if ((c = DD[j] + SS[j] + gi) > midc)
      { midc = c;
        midj = j;
        type = 2;
      }
  if (lc) {
    if ((c = DD[0] + SS[0]) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  } else {
    if ((c = DD[0] + SS[0] + gi) > midc)
      { midc = c;
        midj = 0;
        type = 2;
      }
  }
}
/* Conquer: recursively around midpoint */

  if (midj == 0 || midj == N) {
     if (type == 1)
       { __gps_align(A,     B     ,midi  ,midj  ,tb,-gj,topr,0   ,lc,rc);
         __gps_align(A+midi,B+midj,M-midi,N-midj,-gi,te,0   ,botr,lc,rc);
       }
     else
       { __gps_align(A,       B     ,midi-1  ,midj  ,tb,0 ,topr,0   ,lc,rc);
         DEL(2);
         __gps_align(A+midi+1,B+midj,M-midi-1,N-midj,0 ,te,0   ,botr,lc,rc);
       }
  } else {
     if (type == 1)
       { __gps_align(A,B,midi,midj,tb,-gj,topr,0,lc,0);
         __gps_align(A+midi,B+midj,M-midi,N-midj,-gi,te,0,botr,0,rc);
       }
     else
       { __gps_align(A,B,midi-1,midj,tb,0,topr,0,lc,0);
         DEL(2);
         __gps_align(A+midi+1,B+midj,M-midi-1,N-midj,0,te,0,botr,0,rc);
       }
  }
  return midc;
}

void 
gps_align (A, B, M, N, gopi, gepi, gopj, gepj, ali)
PROFILECOLUMN *A;
SEQTYPE       *B;
int            M, N;
CALCTYPE       gopi, gepi, gopj, gepj;
ALI           *ali;
{

    int i;
    int *S;
    int t;
    CALCTYPE score, s;
    int lasti, lastj, diff, length;

#ifdef DEBUG
    printf("gps_align called with gopi %5.2f, gepi %5.2f, gopj %5.2f, gepj %5.2f, M %i, N %i\n", gopi, gepi, gopj, gepj, M, N);

    PrintProfile ( A, M );
    PrintSequence( B, N );
#endif

    /* Setup global parameters */
    gi = gopi; hi = gepi; mi = gi+hi;
    gj = gopj; hj = gepj; mj = gj+hj;


    S = Malloc( MAXALI_SIZE * sizeof(int) );
    sapp = S;
    last = 0;
    /* do the alignment */
    ali->score = __gps_align(A, B, M, N,
			    -gi, -gj,
			    1,1,1,1);


    /* do some corrections */
    if ((abs(S[0]) < abs(S[1])) && S[0] != 0) {
	t = S[1];
	S[1] = S[0];
	S[0] = t;
    }
    if ((abs(sapp[0]) < abs(sapp[-1])) && sapp[0] != 0) {
	t = sapp[-1];
	sapp[-1] = sapp[0];
	sapp[0] = t;
    }
    
    /* convert to table formatted alignment */
    ConvertAli(M, N, S, ali);

    /* fill in score-array */
    /* fill in score-array and rescore alignment */
    lasti = ali->align_i[0];
    lastj = ali->align_j[0];
    length= ali->lastindex+1;

    for (i = 0; i <= ali->lastindex; i++) {
	/* score matches */
	s = A[ali->align_i[i]][B[ali->align_j[i]]];
	ali->align_s[i] = s;
	score          += s;
	/* penalize gaps */
	if (lasti != ali->align_i[i]) {
	    diff = ali->align_i[i] - lasti;
	    score -= (gi + hi * diff);
	    length += diff;
	}	    
	if (lastj != ali->align_j[i]) {
	    diff = ali->align_j[i] - lastj;
	    score -= (gj + hj * diff);
	    length += diff;
	}	    
	lasti = ali->align_i[i] + 1;
	lastj = ali->align_j[i] + 1;
    }
    
    ali->score  = score;
    ali->length = length;
    
    free (S);

}
