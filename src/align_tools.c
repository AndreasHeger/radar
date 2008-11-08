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
$Id: align_tools.c,v 1.2 2005/02/17 23:13:32 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#include "align.h"
#include "align_tools.h"
#include "toolbox.h"

#define PREVJ { if (--j < 1) j = N; }

void TraceBack ( stopcell, M, N, trace, ali) 
CELL stopcell;
int M, N;
int  * trace;
ALI  * ali;
{
    int i, j;
    int lalign    = 0;
    int lastalign = 0;
    int O = N + 1;
    int t;

    i = stopcell.i;
    j = stopcell.j;

    t = trace[INDEX(i,j)];
    
    while ( t != 0) {
#ifdef SAVE
	if (lastalign > MAXALI_SIZE) {
	    printf("Error: Maximum alignments %i size exceeded; Shut down computer immediately\n", MAXALI_SIZE);
	    exit(1);
	}
#endif
#ifdef DEBUG
	printf("%4i %4i %4i\n", i, j, t);
#endif
	switch (t) {
	case DELETION  : 
	    lalign ++;
	    if (--j < 1) {
		j = N;
		i--;
	    }
	    break;
	case INSERTION : i--; 
	    lalign ++;
	    break;
	case MATCH     : 
	    ali->align_i[lastalign] = i;
	    ali->align_j[lastalign] = j;
	    lastalign++;
	    lalign++;
	    i--; 
	    PREVJ; 
	    break;
	case ROLL :                 /* roll over, with match */ 
	    ali->align_i[lastalign] = i;
	    ali->align_j[lastalign] = j;
	    lastalign++;
	    lalign++;
	    PREVJ;
	    break;
/* 	default: */
/* 	    lalign += N - t; */
/* 	    j = t;          */
/* 	    break; */
	}
	t = trace[INDEX(i,j)];
    } 
#ifdef DEBUG
    printf("\n");
#endif

    ali->length    = lalign;
    ali->lastindex = lastalign-1; 
}

void ConvertAli (M, N, S, ali) 
int      M, N; 
int     *S;
ALI     *ali;
{ 
    register int   i,  j, op;
    int alilength;
    int lastindex;
    
    alilength = i = j = op = 0;
    lastindex = 0;
    
    while (i < M || j < N) {
	op = *S++;
	if (i == 0 && j == 0 && op != 0) {
	    if (op > 0) 
		j += op;
	    else 
		i -= op;
	} else if (i == M || j == N) {
	    i = M;
	    j = N;
	} else if (op == 0) {
	    i++; j++;
	    ali->align_i[lastindex] = i;
	    ali->align_j[lastindex] = j;
	    lastindex++;
	    alilength++;
	} else if (op > 0) {
	    j+= op;
	    alilength += op;
	} else {
	    i-= op;
	    alilength -= op;
	}
    }
    
    ali->lastindex = lastindex - 1;
    ali->length    = alilength;
}

CALCTYPE DotProduct ( A, B )
CALCTYPE *A;
CALCTYPE *B;
{
    int i;
    CALCTYPE score;
    score = 0;
    for (i = 0; i < PROFILEWIDTH; i++ ) {
	score = score + A[i] * B[i];
    }
    return score;
}

CALCTYPE ProfileScore ( A, AC, B, BC)
CALCTYPE  *A;
CALCTYPE *AC;
CALCTYPE  *B;
CALCTYPE *BC;
{
    int i;
    CALCTYPE score;
    score = 0;
    for (i = 0; i < PROFILEWIDTH; i++ ) {
	score += AC[i] * B[i] + BC[i] * A[i];
    }
    return score;
}



int WasGapHere( trace, i, j, M, N )
int  * trace;
int i, j; 
int M, N;
{
    int O = N + 1;
    int t;
    int gapfound = 0;
    int compj    = j;

    j--;
    if (j < 1) {
	j = N;             /* wrap over deletions */
	i--;
    }

    t = trace[INDEX(i,j)];
    while ( (t != 0) && (!gapfound) ) {
	switch (t) {
	case DELETION  : 
	    if (j == compj) {
		gapfound = 1;
	    }
	    j--;
	    if (j < 1) {
		j = N;             /* wrap over deletions */
		i--;
	    }
	    break;
	case INSERTION : i--; 
	    break;
	case MATCH     : 
	    i--; 
	    j--; 
	    if (j < 1) {
		j = N;             /* wrap over deletions */
	    }
	    break;
	case ROLL :                 /* roll over, with match */ 
	    j--;
	    if (j < 1) {
		j = N;
	    }
	    break;
	}
	t = trace[INDEX(i,j)];
    } 

    return gapfound;
}


/*-------------------------------*/
static double q[NCOMPONENTS] = { 0.182962,0.057607,0.089823,0.079297,0.083183,0.091122,0.115962,0.06604,0.234006 } ;

static double a[NCOMPONENTS][PROFILEWIDTH] = {
      { 0.270671, 0.039848, 0.017576, 0.016415, 0.014268, 0.131916, 0.012391, 0.022599, 0.020358, 0.030727, 0.015315, 0.048298, 0.053803, 0.020662, 0.023612, 0.216147, 0.147226, 0.065438, 0.003758, 0.009621 },
      { 0.021465, 0.0103, 0.011741, 0.010883, 0.385651, 0.016416, 0.076196, 0.035329, 0.013921, 0.093517, 0.022034, 0.028593, 0.013086, 0.023011, 0.018866, 0.029156, 0.018153, 0.0361, 0.07177, 0.419641},
      {0.561459, 0.045448, 0.438366, 0.764167, 0.087364, 0.259114, 0.21494, 0.145928, 0.762204, 0.24732, 0.118662, 0.441564, 0.174822, 0.53084, 0.465529, 0.583402, 0.445586, 0.22705, 0.02951, 0.12109 },
      { 0.070143, 0.01114, 0.019479, 0.094657, 0.013162, 0.048038, 0.077, 0.032939, 0.576639, 0.072293, 0.02824, 0.080372, 0.037661, 0.185037, 0.506783, 0.073732, 0.071587, 0.042532, 0.011254, 0.028723 },
      { 0.041103, 0.014794, 0.00561, 0.010216, 0.153602, 0.007797, 0.007175, 0.299635, 0.010849, 0.999446, 0.210189, 0.006127, 0.013021, 0.019798, 0.014509, 0.012049, 0.035799, 0.180085, 0.012744, 0.026466 },
      { 0.115607, 0.037381, 0.012414, 0.018179, 0.051778, 0.017255, 0.004911, 0.796882, 0.017074, 0.285858, 0.075811, 0.014548, 0.015092, 0.011382, 0.012696, 0.027535, 0.088333, 0.94434, 0.004373, 0.016741 },
      { 0.093461, 0.004737, 0.387252, 0.347841, 0.010822, 0.105877, 0.049776, 0.014963, 0.094276, 0.027761, 0.01004, 0.187869, 0.050018, 0.110039, 0.038668, 0.119471, 0.065802, 0.02543, 0.003215, 0.018742 },
      { 0.452171, 0.114613, 0.06246, 0.115702, 0.284246, 0.140204, 0.100358, 0.55023, 0.143995, 0.700649, 0.27658, 0.118569, 0.09747, 0.126673, 0.143634, 0.278983, 0.358482, 0.66175, 0.061533, 0.199373 },
      { 0.005193, 0.004039, 0.006722, 0.006121, 0.003468, 0.016931, 0.003647, 0.002184, 0.005019, 0.00599, 0.001473, 0.004158, 0.009055, 0.00363, 0.006583, 0.003172, 0.00369, 0.002967, 0.002772, 0.002686 }
     };

static double bg[PROFILEWIDTH];

static double wa[PROFILEWIDTH];
                         /*                 a  b c d e f  g h i j k  l  m  n o p  q r s  t  u v  w  x  y  z  */
static unsigned char *encodealphabet    = "ACDEFGHIKLMNPQRSTVWYX-"; 
/* 20 is maskcode */
/* 21 is gapcode */ 

/* coding large characters alphabetically into AA 0 to 19, all others (for example lower case characters)
   are set to GAPCODE */

static SEQTYPE decode[131] = { MASKCODE, 	                                                                       /* 0 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 10 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 20 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 30 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 40 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 50 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 60 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE,       0, MASKCODE,        1,         2,        3,        4, /* 70 */
		            5,        6,        7, MASKCODE,        8,        9,       10,       11, MASKCODE,       12, /* 80 */
		           13,       14,       15,       16, MASKCODE,       17,       18, MASKCODE,       19, MASKCODE, /* 90 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 100 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 110 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, /* 120 */
		     MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE, MASKCODE  /* 130 */
		     };

char Encode( x )
    int x;
{
    if (x <= GAPCODE) 
	return encodealphabet[x];
    else 
	return UNKNOWN_CODE;
}

SEQTYPE * GetDecode() {
    return decode;
}

SEQTYPE Decode( x ) 
    int x;
{ return decode[x]; }

void ProfileInitialize () 
{
    double x[PROFILEWIDTH];
    double total;
    int i,k;
    
    /* calculate null model: uniform aa frequencies -> bg[]*/
    for (i = 0; i < PROFILEWIDTH; i ++) {
	x[i] = 0;
    }
    for (k = 0; k < NCOMPONENTS; k++) {
	for (i = 0; i < PROFILEWIDTH; i ++) {
	    x[i] = x[i] + q[k] * a[k][i];
	}
    }
    
    /* normalize null modell */
    total = 0;
    for (i = 0; i < PROFILEWIDTH; i++) {
	total = total + x[i];
    }
    
    if (total == 0) {
	total = 1;
    }
    
    for (i = 0; i < PROFILEWIDTH; i++) {
	bg[i] = x[i] / total;
    }
    /* calculate |aj| -> wa[] */
    for (k = 0; k < NCOMPONENTS; k++) {
	total = 0;
	for (i = 0; i < PROFILEWIDTH; i ++) {
	    total = total + a[k][i];
	}
	wa[k] = total;
    }

    
}


void RegularizeColumn( count, pc) 
COUNTTYPE  *count;
CALCTYPE   *pc;
{
    int i, k;
    COUNTTYPE ntotal;
    double s;
    double sum1, p;
    double u[NCOMPONENTS];

    ntotal = 0;
    /* calculate number of observations */
    for (i = 0; i < PROFILEWIDTH; i ++) {
	ntotal += count[i];
    }

    /* precalculate some terms */
    s = 0;
    for (k = 0; k < NCOMPONENTS; k ++) {
	sum1 = 0;
	for (i = 0; i < PROFILEWIDTH; i ++) {
	    sum1 += (double)count[i] * a[k][i]; /* count_i * alpha_ki */
	}
	u[k] = q[k] * sum1;		/* u = mixture_coeff * sum over aa of count_i * alpha_i */
	s = s + u[k];			/* s = sum over all Dirichlet distribution(n) */
    }


    /* if there were no observations (s == 0) */
    if (s == 0) {
	s = NCOMPONENTS;
	for (k = 0; k < NCOMPONENTS; k++) {
	    u[k] = 1;			/* each component is equally likely */
	}
    }

    /* calculate posterior probability for each aa in Profilecolumn */
    for (i = 0; i < PROFILEWIDTH; i ++) {
	p = 0;
	for (k = 0; k < NCOMPONENTS; k ++) {
	    p = p + u[k] / s * ( (double)count[i] + a[k][i]) / ((double)ntotal + wa[k] );
	}

	/* log odds score */
	if (p <= 0) {
	    p = 0;
	} else {
	    p = log(p / bg[i]);
	}
	
	/* save*/

	pc[i] = (CALCTYPE)p;
    }
}    

void RegularizeColumnProbability( count, pc) 
COUNTTYPE  *count;
CALCTYPE   *pc;
{
    int i, k;
    COUNTTYPE ntotal;
    double s;
    double sum1, p;
    double u[NCOMPONENTS];

    ntotal = 0;
    /* calculate number of observations */
    for (i = 0; i < PROFILEWIDTH; i ++) {
	ntotal += count[i];
    }

    /* retreat, if not enough observations and set column to zeros */
    if (ntotal == 0) {
	for (i = 0; i < PROFILEWIDTH; i ++) 
	    pc[i] = 0;
	return;
    }

    /* precalculate some terms */
    s = 0;
    for (k = 0; k < NCOMPONENTS; k ++) {
	sum1 = 0;
	for (i = 0; i < PROFILEWIDTH; i ++) {
	    sum1 += (double)count[i] * a[k][i]; /* count_i * alpha_ki */
	}
	u[k] = q[k] * sum1;		/* u = mixture_coeff * sum over aa of count_i * alpha_i */
	s = s + u[k];			/* s = sum over all Dirichlet distribution(n) */
    }


    /* if there were no observations (s == 0) */
    if (s == 0) {
	s = NCOMPONENTS;
	for (k = 0; k < NCOMPONENTS; k++) {
	    u[k] = 1;			/* each component is equally likely */
	}
    }

    /* calculate posterior probability for each aa in Profilecolumn */
    for (i = 0; i < PROFILEWIDTH; i ++) {
	p = 0;
	for (k = 0; k < NCOMPONENTS; k ++) {
	    p = p + u[k] / s * ( (double)count[i] + a[k][i]) / ((double)ntotal + wa[k] );
	}
	
	/* save (i.e. skip the log-odds-step) */
	pc[i] = (CALCTYPE)p;
    }

#ifdef SAVE
    /* well, they don't always sum up to 1. Why? probably the neglection of the Gamma-Function-Terms */
    /* normalize, so we do have probabilities */
    p = 0;
    for (i = 0; i < PROFILEWIDTH; i++) 
	p += pc[i];
    for (i = 0; i < PROFILEWIDTH; i++) 
	pc[i] = pc[i] / p;
#endif

}    


void Counts2Frequencies ( count, M, frequency )
COUNTCOLUMN     *count;
int M;
FREQUENCYCOLUMN *frequency;
{
    int i, j;
    int ntotal;

    for ( j = 1; j <= M; j++) {
	ntotal = 0;
	/* calculate number of observations */
 	for (i = 0; i < PROFILEWIDTH; i ++) {
	    ntotal += count[j][i];
	}
	/* copy into frequency column */
	for (i = 0; i < PROFILEWIDTH; i ++) {
	    frequency[j][i] = (double)count[j][i] / (double)ntotal;
	}
    }
}


void QuickSortIndices (from, to, x, y)
int from, to;
int *x, *y;
{
    int lastsmall, t, c1, i;
    int m; /* value of median */

/*     printf(" from %i, to %i \n", from, to); */
/*     for (i = 0; i < 17; i++) printf("%i - ", y[x[i]]); */
/*     printf("\n"); */
    if (from < to ) {
	c1 = (from + to) / 2;         
	t = x[from]; x[from] = x[c1]; x[c1] = t;
	m = y[x[from]]; /* choose median value to compare to */
	lastsmall = from;
	for (i = from + 1; i <= to; i++) {
	    if ( y[x[i]] < m) { /* swap lastsmall and i */
		lastsmall++;
		t = x[lastsmall]; x[lastsmall] = x[i]; x[i] = t;
	    }
	}
	
	t = x[from]; x[from] = x[lastsmall]; x[lastsmall] = t; 
	
	m = lastsmall;
	QuickSortIndices( from, m, x, y);
	QuickSortIndices( m+1, to, x , y);
    }
}

/* sort dots by row and then by column */
void SortDots( row, col, score, ndots ) 
    int *row, *col;
CALCTYPE *score;
int ndots;
{

    int       i,x, from, to;
    int      *temp1;
    CALCTYPE *temp2;
    int      *rowindices = (int*)Malloc( ndots * sizeof(int));

    /* sort indices on row */

    for (i = 0; i < ndots; i++) {
	rowindices[i] = i;
    }

    QuickSortIndices( 0, ndots - 1, rowindices, row);

    /* sort indices on column inside row */
    from = 0;
    while ( from < ndots ) {
	x    = row[rowindices[from]];
	to   = from + 1;
	while ( (to < ndots) && (x == row[rowindices[to]]) ) { to++; } /* find end of row */
	QuickSortIndices( from, to - 1 , rowindices, col ); /* and sort per column */
	from = to;
    }
    
    /* now reorder; takes a lot of space  and time, improve later*/
    temp1 = Malloc( ndots * sizeof( int ));
    
    for (i = 0; i < ndots; i++) {
	temp1[i] = row[rowindices[i]];
    }
    memcpy( row, temp1, ndots * sizeof (int));

    for (i = 0; i < ndots; i++) {
	temp1[i] = col[rowindices[i]];
    }
    memcpy( col, temp1, ndots * sizeof (int));

    temp2 = Malloc( ndots * sizeof( CALCTYPE ));
    for (i = 0; i < ndots; i++) {
	temp2[i] = score[rowindices[i]];
    }
    memcpy( score, temp2, ndots * sizeof (CALCTYPE));

    free (rowindices);
    free (temp1);
    free (temp2);
}

/* Example:
row = 1 , 2, 3, 4, 5, 6, 3, 4, 5, 6
col = 10,11,12,13,14,15,20,21,22,23
is converted to x1
x1  = 1,2,3,7,4,8,5,9,6,10
returns:
firstdot = 1
nextrow  = 2,3,7,8,9,10,4,5,6,-1
*/

/* build array for fast access to dots */
/* each element contains index in row for sequence position i */
/* note: row has to be already sorted in ascending order !! */
/* memory for dotp has to be allocated before */ 
int * CalculateDotp ( row, ndots, maxrow, dotp )
int * row;
int ndots;
int maxrow;
int * dotp;
{
    int firstdot = 0;
    int thisdot;
    int i;
    
    for (i = 0; i <= maxrow; i++) { dotp[i] = -1; } /* clear */
    
    for (i = 0; i < ndots; i++) {
	thisdot = row[i];
        if (thisdot != row[firstdot]) firstdot = i;
	if( thisdot > maxrow )        break; /* exit, if beyond limit */
	dotp[thisdot] = firstdot;
    }

#ifdef DEBUG 
    /*    printf("CalculateDotp\n");
    for (i = 1; i <= maxrow; i++) {
	printf("%i dotp[i] %i\n", i, dotp[i]);
	}*/
#endif
    return dotp;
}

/*--------------------------------------------------------------------*/
/* retrieves the index for a dot given by i and j. If the dot does not
   exist, then it returns -1 */
int GetDotIndex( r , c , row, col, dotp)
    int r,c;
int *row;
int *col;
int *dotp;
{                     
    int x     = dotp[r];
    int found = 0;

    if (x < 0) return -1;

    while (row[x] == r ) {
        if (col[x] == c) {
	    found = 1;
            break;
        }
        x++;
    }
    
    if (found)
        return x;
    else
        return -1;
}

void TraceBackDotalignment (
	       row, col, 
	       score, ndots,
	       roll, 
	       minrow, maxrow, 
	       mincol, maxcol, 
	       lastdot, 
	       trace,
	       ali
	       )
int     *row;
int     *col;
double  *score;
int      ndots;
int      roll;                  
int      minrow, maxrow; 
int      mincol, maxcol;
int      lastdot;
int     *trace;
ALI     *ali;
{
    int last_i, last_j, alilength;
    int jleft, idot, i, ires, jres;

    /* trace back ----------------- */
#ifdef DEBUG
    printf("TraceBackDotalignment starts\n");
#endif
    /* print traceback */
    idot  = lastdot; 
    jleft = maxrow;
    i = 0;
    alilength = 0; last_i = col[idot]; last_j = row[idot];
    
    while ( idot >= 0) { 

#ifdef DEBUG
	printf("-->idot %i; col[idot] %i; row[idot] %i; trace[idot] %i\n",
	       idot,col[idot],row[idot],trace[idot]) ; 
#endif
	jres = row[idot]; 
	ires = col[idot];
	ali->align_s[i] = score[idot];

	if (jres < 1) continue; 
	if (ires < 1) continue; /*  hacky bug fix */
	if ( jres > jleft) break; 
	jleft = jres; /* just in case */
	if (roll) {                   /* when rolling, row and columns were exchanged, I change this back */
	    ali->align_i[i] = ires;
	    ali->align_j[i] = jres;
	    if (last_i > ires) 
		alilength += last_i - ires - 1;
	    else if (last_i != ires)
		alilength += last_i - mincol + maxcol - ires;
	    if (last_j > jres) 
		alilength += last_j - jres - 1;
	    else if (last_j != jres)
		alilength += last_j - minrow + maxrow - jres;
	} else {	/* map sbjct to repeat */
	    ali->align_i[i] = ires;
	    ali->align_j[i] = jres;
	    if (last_i != ires) alilength += last_i - ires - 1;
	    if (last_j != jres) alilength += last_j - jres - 1;
	}		/* map query to sbjct */
	last_i = ires;
	last_j = jres;
	i++;
	idot = trace[idot];
    }
    ali->length    = alilength + i;
    ali->lastindex = --i;
    
}

int GetFirstDot( jres, dotp ) 
    int jres;
    int *dotp;
{
    return dotp[jres];
}
