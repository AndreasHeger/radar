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
$Id: dotalign.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/
#include "align.h"
#include "align_tools.h"
#include "toolbox.h"
/* quicksort x, which contains the indices for y, which
   contains the values to be sorted 
*/

/* this array starts at 1, since it maps sequence to dotplot */
static int *dotp;

/* assumptions: row, col and score start at 0
   ndots: number of dots = lastelement(row) +1
   in row and column: aa are labeled starting at 1
*/

/* 
   hdot: dot of last iteration
   idot: dot of current iteration
   xdot: dot for local iteration
*/

CALCTYPE __DotAlign ( row, col, 
		      score, ndots,
		      roll, 
		      minrow, maxrow, 
		      mincol, maxcol, 
		      jopen, jelon, 
		      iopen, ielon, 
		      lastdot, 
		      trace )
int     *row;
int     *col;
double  *score;
int      ndots;
int      roll;                  
int      minrow, maxrow; 
int      mincol, maxcol;
CALCTYPE jopen, jelon;
CALCTYPE iopen, ielon;
int     *lastdot;
int     *trace;

{
    
    int left;
    int nheu, nadg;
  
    int bestdot, globdot;
    CALCTYPE globbest, best;
    CALCTYPE sa,sb,sc,sd,se,sf,s;

    int i;

    int acol, xcol, fcol;
    int adot, bdot, cdot, ddot, edot, fdot, hdot, idot, xdot;
    int ires, jres;
    
    int bestpercolstackptr; /* points to next free element in bestpercolstack */
  
  /* Malloc for bestpercol and topdot */
    int *bestpercolstack = (int*)Malloc( (maxcol + 1) * sizeof (int));
    /* 7.1.00, update other, per column has to have space for complete row */
    int *bestpercol      = (int*)Malloc( (maxrow + 1) * sizeof (int)); 
    int *topdot          = (int*)Malloc( (maxrow + 1) * sizeof (int));
    CALCTYPE *m          = (CALCTYPE*)Malloc( ndots * sizeof(CALCTYPE));

    /* initialisation */
    adot    =-1; bdot    =-1; cdot =-1; ddot =-1; edot =-1; fdot=-1; hdot=-1;
    globdot =-1; globbest =0; best = 0; nadg = 0; nheu = 0;
    bestpercolstackptr = STACKEMPTY;
    
    for (ires = mincol; ires <= maxcol; ires++) { bestpercol[ires] = -1; }
    for (jres = minrow; jres <= maxrow; jres++) { topdot[jres]     = -1; }
    for (i = 0; i < ndots; i ++) { trace[i]     = -1; m[i] = 0; }
    acol = mincol;
    idot = 0;
    
  /* iterate through nextrow */
    for ( idot = dotp[minrow]; idot < ndots; idot++ ) { /* iterate through nextrow starting at firstrow */
	if (idot < 0) continue;
	jres = row[idot];                           /* jres = row */
	ires = col[idot];                           /* ires = col, wrap around col */

	if (jres > maxrow) break;                   /* checking boundaries */
	if (ires < mincol) continue;                /* can be removed, if only dots in boundaries are supplied */
	if (ires > maxcol) continue;

#ifdef DEBUG
	printf("--------------------------------------------\n");
	printf("idot = %i, jres = %i, ires = %i, score = %5.2f\n", idot, jres, ires, score[idot]);
#endif
	/* calculate top row */
	if ( (hdot < 0) ||           /* enter first time */
	     (jres > row[hdot]) ) {  /* skip, if not in the same row as last time*/
	    
	    topdot[jres] = idot;
	    while( bestpercolstackptr > STACKEMPTY ) {
		xdot = bestpercolstack[--bestpercolstackptr];
		if ( row[xdot] >= jres ) break;                 /* stop, if entering current row */
		xcol = col[xdot];
		if (xdot < 0 )            continue;             /* safety check */
		if (bestpercol[xcol] < 0) 
		    bestpercol[xcol] = xdot;
		else if (m[xdot] > m[bestpercol[xcol]] ) 
		    bestpercol[xcol] = xdot;
	    }
	    /* init all */
	    acol = mincol; fcol = ires + 1; 
	    adot = -1; bdot = -1; cdot = -1; edot = -1; fdot = -1; 
	}
	/* update c = maximum dot along last row */
	xdot = dotp[jres - 1];
	sc = 0; 
	while ( (xdot > -1)  && 
		(row[xdot] == jres - 1 ) &&  /* stop, if dot in previous row any more*/
		(col[xdot]  < ires - 1 )     /* end, if direct contact to new dot*/
		) { 

	    s = m[xdot] + iopen + ielon * (ires - col[xdot] - 1);
	    
	    if (s > sc) { 
		cdot = xdot; 
		s    = sc; 
	    }
	    xdot++;
	}
	
	/* compute d = match adjacent dot in previous row and column */
	if (roll) { 
	    if ( ires == mincol) { 
		i = maxcol; 
	    } else { 
		i = ires - 1; } 
	} else { 
	    if (ires > mincol) { 
		i = ires - 1; 
	    } else { 
		i = -1; 
	    }
	}
	
        /* ddot -> dot in previous column, previous row */
	ddot = GetDotIndex( jres - 1, i, row, col, dotp ); /* look up index in row, col, score for dot; -1 if not found */
	/* update a = */
	
	if ( ddot < 0 ) { /* only update a if d unoccupied */

	    if ( adot > -1 && xdot > -1 ) { /* Hmm. previous match? */ 
		sa = m[adot] + 
		    jopen + jelon * (jres - row[adot] - 1) +
		    iopen + ielon * (ires - col[xdot] - 1); 
	    } else { 
		sa=0; 
	    }
	    
	    /* ##>> searcharea! */
	    left = 0;
	    for (i = acol; i <= ires-2; i++) {
		xdot = bestpercol[i]; 
		if ( xdot < 0 ) continue;
		s = m[xdot] + 
		    jopen + jelon * (jres - row[xdot] - 1) + 
		    iopen + ielon * (ires - col[xdot] - 1);
		if (s > sa) { 
		    sa=s; 
		    adot=xdot; 
		}
		if (row[xdot] > left) { 
		    left = row[xdot]; 
		}
		if(left == jres-2) break;
	    }
	    acol = ires;
	    nadg++;
	} else { 
	    adot = -1; 
	    nheu++; 
	}
	/* update bdot = max scoring dot in previous column */

	if ( roll ) { 
	    if ( ires==mincol ) { 
		bdot = bestpercol[maxcol]; 
	    } else { 
		bdot = bestpercol[ires-1]; 
	    } 
	} else { 
	    if( ires > mincol) { 
		bdot = bestpercol[ires-1]; 
	    } else { 
		bdot = -1; 
	    } 
	}
	if (bdot == ddot) bdot = -1; /* dot has to be in different row */
	
	/* rollover */
	if(roll) {
	    /* update f */
	    if ( fdot > -1 ) { 
		sf = m[fdot] + 
		    jopen + jelon * ( jres   - row[fdot] - 1 ) + 
		    iopen + ielon * ( maxcol - col[fdot] + ires - mincol); 
		if( col[fdot] <= ires) 
		    fcol = ires + 1; 
	    } else { 
		sf=0; 
	    }
	    /*##>> searcharea! */
	    left = 0;

	    for ( i = fcol; i <= maxcol; i++) { /*# only executed for topdot or f<=ires */
		if ( i < 0 ) continue;          /* safety */
		xdot = bestpercol[i]; 
		if ( xdot < 0 ) continue;
		s = m[xdot] + 
		    jopen + jelon * (jres   - row[xdot] - 1) + 
		    iopen + ielon * (maxcol - col[xdot] + ires - mincol);
		if ( s > sf ) { 
		    sf   = s; 
		    fdot = xdot;
		}

		if (row[xdot] > left) {
		    left = row[xdot];
		}

		if ( left == jres-2 ) break;
	    }

	    fcol = maxcol + 1;
	    /* update e = */
	    if ( (edot > -1 && col[edot] <= ires) || (edot < 0) ) { 
		if ( jres > mincol ) {
		    xdot = topdot[jres-1];
		} else {
		    xdot =- 1;
		}
		se   =  0; 
		edot = -1;
		if (xdot >= 0) {
		    for ( xdot = 0; xdot < ndots; xdot++ ) {
			if (row[xdot] != jres - 1 ) break;	/* since sorted by row first. check, if we leave the row */
			if (col[xdot] < mincol) continue;	/* since sorted by column inside a column */
			if (col[xdot] > maxcol) break;		
			if (col[xdot] >  ires) {
			    s = m[xdot] + iopen + ielon * (ires - mincol + maxcol - col[xdot]);
			    if (s>se) {
				se=s; 
				edot=xdot;
			    }
			}
		    }
		}
	    } 
	} /* endof if(roll)		*/
	
	/*  select best of d|a|b|c */
	best = 0; bestdot = -1; 
	sa = sb = sc = sd = se = sf = 0;
	if ( ddot >= 0) {	    
	    sd = m[ddot];
	    if( sd > best ) { best = sd; bestdot = ddot; }
	} 
	if ( adot >= 0 ) {
	    sa = m[adot] + jopen + jelon * (jres - row[adot] - 1) + iopen + ielon * (ires - col[adot] - 1); 
	    if( sa > best) { best = sa; bestdot = adot; }
	}
	if ( bdot >= 0 ) {
	    sb = m[bdot] + jopen + jelon * (jres - row[bdot] - 1); 
	    if( sb > best) { best = sb; bestdot = bdot; }
	}
	if ( cdot >= 0) {
	    sc = m[cdot] + iopen + ielon * (ires - col[cdot] - 1); 
	    if( sc > best) { best = sc; bestdot = cdot; }
	}
	if (roll) {
	    if (edot >= 0) {
		se = m[edot] + iopen + ielon * (ires - col[edot] + maxcol - mincol); 
		if( se > best ) { best = se; bestdot = edot; }
	    }

	    if (fdot >= 0) {
		sf = m[fdot] + jopen + jelon * (jres - row[fdot]) + iopen + ielon * (ires - col[fdot] + maxcol - mincol); 
		if( sf > best ) { best = sf; bestdot = fdot; }
	    }
	}
	/* record traceback */
	best += score[idot]; 
	
	if (best < 0) { /* local alignment, reset to zero or start new trace with single match */
	    best    = 0 ; 
	    bestdot = -1;
	}				
	m[idot]     = best;
	trace[idot] = bestdot; 
	
	if ( best > globbest) { /* save best dot */
	    globbest = best; 
	    globdot  = idot;
	}
	
#ifdef DEBUG
	printf("idot %5i; trace[idot] %5i; m[idot] %5.2f; best %5.2f\n",
	       idot, trace[idot], m[idot], best);
	printf("adot %5i; bdot %5i; cdot %5i; ddot %5i; edot %5i; fdot %5i\n",
	       adot,bdot,cdot,ddot,edot,fdot );
	printf("sa   %5.2f; sb   %5.2f; sc   %5.2f; sd   %5.2f; se   %5.2f; sf   %5.2f\n",
	       sa, sb, sc, sd, se, sf);
	printf("acol %5i; globdot %5i; globbest %5.2f\n",
	       acol, globdot, globbest );
#endif
	/* for bestpercol */
	
	if ( (bestpercol[ires] == -1) ) 
	    bestpercolstack[bestpercolstackptr++] = idot;
	else if (best > m[bestpercol[ires]]) 
	    bestpercolstack[bestpercolstackptr++] = idot;

	if (bestpercolstackptr > maxrow) 
	    printf ("bestpercolstackptr %i idot %i %i %i\n", bestpercolstackptr, idot, row[idot], col[idot]);

	/*  h=i */
	hdot = idot;
	/* #if((nadg+nheu) =~ /00/) { warn "interim nheu=nheu nadg=nadg\n"; } */
    }
    
#ifdef DEBUG  
    printf("total nheu= %i nadg= %i\n", nheu, nadg);
#endif

    *lastdot = globdot;
    
    free( bestpercolstack );
    free( bestpercol );        
    free( topdot );                
    free( m );        
    
    return ( globbest );
}

void dotalign ( row, col, score, 
		ndots,
		roll, 
		minrow, maxrow, /* changed: 4.1.00, */
		mincol, maxcol, /* changed: 4.1.00 */
		iopen, ielon, jopen, jelon, 
		ali )
int     *col;
int     *row;
double  *score;
int      ndots;
int      roll;                  
int      minrow, maxrow, mincol, maxcol; 
CALCTYPE iopen, ielon, jopen, jelon;
ALI     *ali;

{ 

    int lastdot;
    int *trace = (int*)Malloc( ndots * sizeof(int));
    
#ifdef DEBUG
    printf("DotAlign called with mincol %i, maxcol %i , minrow %i, maxrow %i, roll %i, ndots %i\n",
	   mincol, maxcol, minrow, maxrow, roll, ndots) ;
    printf("gap parameters: iopen %5.2f, ielon %5.2f, jopen %5.2f, jelon %5.2f\n",
	   iopen, ielon, jopen, jelon);
    PrintDots(row, col, score, ndots); 
#endif

    SortDots( row, col, score, ndots);

#ifdef DEBUG
    printf("After sorting:\n");
    /*    PrintDots(row, col, score, ndots); */
#endif
  
    dotp = (int*)Malloc( (maxrow + 2) * sizeof(int) ); /* allocate memory */
    CalculateDotp( row, ndots, maxrow, dotp ); 

    ali->score = __DotAlign ( 
			     row, col, 
			     score, ndots,
			     roll, 
			     minrow, maxrow, 
			     mincol, maxcol, 
			     jopen, jelon, 
			     iopen, ielon, 
			     &lastdot, 
			     trace 
			     );
    
    TraceBackDotalignment (
	       row, col, 
	       score, ndots,
	       roll, 
	       minrow, maxrow, 
	       mincol, maxcol, 
	       lastdot, 
	       trace,
	       ali
	       );
    
    free ( trace);
    free ( dotp );
}
