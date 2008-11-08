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
$Id: register.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository

*/
#include "radar.h"

static CALCTYPE rg_gapopen_i = -1;   /* gap penalty for insertions in profile */
static CALCTYPE rg_gapelon_i = -0.25; 
static CALCTYPE rg_gapopen_j = -4;   /* gap penalty for insertions in sequence */
static CALCTYPE rg_gapelon_j = -1; 
static CALCTYPE rg_gapopen = -12;	/* gap penalty for shiftali */
static CALCTYPE rg_gapelon = -2; 

static double exploreratio = 0.40; /* how many residues to explore around repeat ends for cutting, 20 % of repeatlength on either side */
static int compfactor      = 5;    /* when to assume, that ends are noise, the larger, the more resctrictive */   
static int splitdistance        = 50;   /* residue-separation, above which two repeats are regarded as distinct */

extern VERBOSITY verbose;

/*-------------------------------------------------------------*/
SEQUENCEWINDOW GetFirstRepeat ( ali )
    ARRAYALI * ali;
{
    
    int to, from;
    SEQUENCEWINDOW w;
    int val;
    int lsequence = ali->length;

    from = 1;
    while( ali->ali[from] < 1 && from <= lsequence ) { from++; } 
    val = ali->ali[from];                      
    to = from + 1;
    
    while( (to <= lsequence) &&  
	   ((ali->ali[to] < 1   ) ||			/* skip gaps */
	    (ali->ali[to] > val))			/* find before next repeat */
	   ) { 
	if (ali->ali[to] > 0) {
	    val = ali->ali[to];
	}
	to++; 
    }     
    to--;
    
    while( ali->ali[to] < 1) { to--; }			/* walk back gaps */
    
    w.from = from;
    w.to   = to;
    return (w);
}

/*-------------------------------------------------------------*/
SEQUENCEWINDOW GetLastRepeat ( ali )
    ARRAYALI *ali;
{
    int to, from;
    SEQUENCEWINDOW w;
    int val;
    int lsequence = ali->length;

    to = lsequence;
    while( ali->ali[to] < 1 && to > 1 ) { to--; } 
    val = ali->ali[to];
    from = to - 1;
    while( ( from > 1) && 
	   ( (ali->ali[from] < 1) ||			/*  skip gaps */
	     (ali->ali[from] < val) )			/*  find before next repeat */
	   )  { 
	if (ali->ali[from] > 0) {
	    val = ali->ali[from];
	}
	from--; 
    }     
    from++;
    
    while ( ali->ali[from] < 1) { from++; }		/* walk back gaps */
  
    w.from = from;
    w.to   = to;
						  
    return (w);
}


/*----------------------------------------------------------
 Title: ShiftAli
 Function: Shifts an alignment, so that residue startcol is
           in the first column and scores the repeats.
           performs no elimination of repeats
           Inserts between repeats are NOT penalized as well as
           incomplete repeats at the end and beginning.
 Returns:  new alignment
           repeatscore = ref to array of scores for each repeat
           sumscore    = sum of scores of all repeats
           nrepeats    = number of repeats
----------------------------------------------------------*/

CALCTYPE ShiftAli ( startcol, mincol, maxcol, 
		    ali, 
		    repeatscores,
		    nrepeats)

int startcol, mincol, maxcol;
ARRAYALI *ali;
CALCTYPE *repeatscores;
int *nrepeats;
{
    int i, j;
    CALCTYPE sumscore, s;
    int lastres, newres, lrepeat, ngaps;
    int lnrepeats;
    int lali = ali->length;
    
    int * map = (int*)Malloc((maxcol + 1) * sizeof(int));

    /* build map */
    j = 0;
    for (i = 0; i <= maxcol; i++) map[i] = 0;
    for (i = startcol; i <= maxcol;   i++) 
	if (ali->ali[i] != 0) 
	    map[i]=++j; 
    for (i = mincol;   i <  startcol; i++) 
	if (ali->ali[i] != 0)
	    map[i]=++j; 
    
    /* go to beginning of first repeat */
    i = 1; while (ali->ali[i] == 0) i++;

    /* intialize */
    lnrepeats = 0;
    lrepeat  = maxcol - mincol + 1;
    lastres  = 0;
    newres   = map[ali->ali[i]];

    sumscore = 0;

    /* compute gap penalty for first repeat */
    if (newres > 1) 
	s = rg_gapopen + rg_gapelon * (newres - 1);
    else 
	s = 0;
    
    /* start main loop */
    while (i <= lali) {
	
	ngaps = 0;
	while ((i <= lali) && ali->ali[i] == 0) { 	           /* count gaps and advance to next non-gap character */
	    ngaps++;
	    i++;
	}
	if (i > lali) break;

	newres = map[ali->ali[i]];
	

	if (newres < lastres) {                                    /* end previous repeat and start a new repeat */
	    if (lrepeat - lastres > 0) {	                   /* penalize for deletions in profile at the end of the repeat */
		s += rg_gapopen + rg_gapelon * (lrepeat - lastres );
	    }

	    repeatscores[lnrepeats++] = s;                          /* save score */
	    sumscore += s;					    /* and add to total */

	    if (newres > 1) 	                                    /* compute gap penalty for start of repeat */
		s = rg_gapopen + rg_gapelon * (newres - 1);
	    else 
		s = 0;
	} else {                                                    /* else penalize for gaps if continuation of same repeat*/
	    if (ngaps) s += rg_gapopen + ngaps * rg_gapelon;
	}

	s += ali->scores[i];                                        /* add alignment score to repeat score */
	ali->ali[i]  = newres;                                      /* add residue to repeat */
	i++;                                                        /* advance to next residue */
	lastres  = newres;                                        
    }

    if (lrepeat - lastres > 0) 
	s += rg_gapopen + rg_gapelon * (lrepeat - lastres );
    
    repeatscores[lnrepeats++] = s;
    sumscore += s;

    *nrepeats = lnrepeats;
    free( map );
    return sumscore;
}


/*------------------------------------------------------------------------
 ----------------------------------------------------------------------
  Title: GetGapVote
  Function:   calculate internal gap lengths
              go through the alignment. For each position in the alignment
              calculate the number of gaps, that occur before that position.
  Parameters: $from: where to start looking
              $to:   where to stop
              $origali: alignment
  Returns:    Hash
  Example:    deletions in sequence
              1234567890
              xx___xxxxx
              xxx_xxxxxx
              xxxx_xxxxx
              xxxx_xxxxx
  Example:    insertions in sequence
              123  45 67890
              xxxyyxx xxxxx insert 1
              xxx  xx xxxxx no insertion
              xxx  xxyxxxxx insert 2
               
  GapVotes: 1,2,3,4,7, 8, 9, 0 = 0;
            5 = 1,
            6 = 5.
  i.e. a good cut would be before residue position 6.
 new version: calculate penalties and intermediate values
 ----------------------------------------------------------------------
*/

CALCTYPE * GetGapVote (from, to, ali, lrepeat ) 
int from, to;
ARRAYALI *ali;
int lrepeat;
{
    int i;
    int thisres, lastres;
    CALCTYPE *gapvote = (CALCTYPE*)Malloc( (ali->length + 1) * sizeof (CALCTYPE));
    
    int max = 0;
    int nins, xres;
    int gaplen;

    for (i = from; i <= to; i++) 
	if (ali->ali[i] > max) max = ali->ali[i] ;

    for (i = 0; i <= ali->length; i++) gapvote[i] = 0;

    /* Initialisierung */
    lastres = ali->ali[to];
    i = to - 1;
    while (i > from ) {
	nins  = 0;							/* reset number of gaps */
	while ( (ali->ali[i] == 0)  && ( i > from )) { nins++; i--; };	/*iterate through gaps in sequence */
	thisres = ali->ali[i];                                          /*thisres: this aligned residue */
	
	if (nins) 
	    gapvote[lastres] += rg_gapopen_i + rg_gapelon_i * nins;     /*penalize for insertions in sequence */
	/*deal with deletions in the sequence, i.e. insertions in the model-repeat  */
	xres   = lastres - 1;                  
	if (xres < 1) xres = max;		             		/*wrap around  */
	if (xres != thisres) {
	    gaplen = xres - thisres;		        		/*gap-penalty for insertions in profile  */
	    if (gaplen < 0) gaplen += max; 				/*modify gap penalty for wrapping around  */
	    while (xres != thisres) {
		gapvote[xres] += rg_gapopen_j + gaplen * rg_gapelon_j;	/*penalize for deletions in sequence  */
		if (nins) 
		    gapvote[xres] += rg_gapopen_i + rg_gapelon_i * nins;/*penalize for insertions in sequence  */
		xres --; gaplen --;
		if (xres < 1) xres = max;				/*wrap around  */
	    }
	}
	
	/* do not decrease nins, since when there is a cut in an insert, the whole insert vanishes */
	
	lastres = thisres;
	i--;                                                                  /*next residue  */
    }
    
    return (gapvote);
}


/* assume that profile starts at 1 */

CALCTYPE  * GetGapPenalties ( ali, lprofile ) 
ARRAYALI *ali;
int lprofile;
{

    int newcol, lastcol;
    int i,j;
    int ngaps;
    CALCTYPE * profilepenalties = Malloc( (lprofile + 1) * sizeof(CALCTYPE ));

    for (i = 0; i <= lprofile; i++) profilepenalties[i] = 0;

    /* advance to first aligned residue */
    i = 1; while ( i <= ali->length && ali->ali[i] == 0) i++;

    lastcol = lprofile;
    newcol = ali->ali[i];
    while ( i <= ali->length) {

	/* calculate for insertions in profile */
	if (newcol < lastcol) { 	/* wrapping around case, insertions in profile */
	    if ( lprofile - lastcol > 0) 
		for (j = lastcol + 1; j <= lprofile; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( lprofile - lastcol) ;
	    if (newcol > 1)
		for (j = 1; j < newcol; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( newcol - 1) ;
	} else {                           
	    if (newcol > lastcol + 1 ) { /* insertion in profile */
		for (j = lastcol + 1; j < newcol; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( newcol - lastcol - 1 ) ;
	    }
	}
	
	lastcol = newcol;
	i++;
	ngaps = 0; while ( i <= ali->length && ali->ali[i] == 0) { i++; ngaps ++; }
	if (i > ali->length)
	    break;

	newcol = ali->ali[i];

	if (ngaps != 0) {
	    profilepenalties[newcol] += rg_gapopen_j + rg_gapelon_j * ( ngaps );
	}
    }
    
    if (verbose > LL3) 
	for (i = 1; i <= lprofile; i++) printf("%i %5.2f\n", i, profilepenalties[i]);

    return (profilepenalties);
}

/*-------------------------------------------------------------------
  after transformation, this means cutting in non-overlapping region, which
  extends from last.to to first.to
  returns: the best profile-column in ali
*/

int CutInD ( window, ali, lrepeat ) 
    SEQUENCEWINDOW window;
ARRAYALI *ali;
int lrepeat;
{
    CALCTYPE bestscore = 0;
    int bestcol        = 0; /* default: no position above 0 found */
    CALCTYPE voted;
    int j = 0;
    
    CALCTYPE * gappenalties = GetGapPenalties( ali, lrepeat );
    
    if (verbose > LL2) 
	printf("Cutting in D called; analysing profile from %i to %i\ni\tgapvote\n", window.from, window.to);
    
    for (j = window.from; j <= window.to; j++) {
	
	voted = -gappenalties[j];          /* the more negative, the better, since the more gaps get annihilated */
	if (voted > bestscore) {           /* >, so the first of equivalent position is chosen */
	    bestscore = voted;
	    bestcol   = j;
	}
	if (verbose > LL2) 
	    printf ("%i\t%5.2f\n", j, voted);
    }
    
    free ( gappenalties );
    
    if (verbose > LL2) 
	printf("Selected: bestcol %i with score %5.2f\n", bestcol, bestscore);
    
    return (bestcol);
}

/*-------------------------------------------------------------------
  after transformation, this means cutting in overlapping region, which
  extends from 1 to last.to
  returns: the best column in ali
------------------------------------------------------------------------*/

int CutInE (first, last, ali, lrepeat) 
SEQUENCEWINDOW first;
SEQUENCEWINDOW last;
ARRAYALI *ali;
int lrepeat;
{

    /* score maximum columns against profile # note: approximation */
    /* move forward the boundary, and remember maximum score */
    int i = first.to;
    int j = last.to;
    CALCTYPE s = 0;
    int bestcol          = first.to + 1;	/* start: non-overlapping alignment; careful, might be gap */
    CALCTYPE bestscore   = 0;
    CALCTYPE lastgapcost = 0;           
    CALCTYPE gapcost;

    CALCTYPE *gapvote = (CALCTYPE*)GetGapPenalties( ali, lrepeat);
    if (verbose > LL2) 
	printf("Case VPG-1:\n scorei\tscorej\ti\tj\tali_i\tali_j\ts\tgapcost\tlastgapcost\n");

    while (i >= first.from) { 

	/* calculate difference scores for residues */
	/* move through gaps */
	while (ali->ali[i] == 0 && i >= first.from ) i--;
	while (ali->ali[j] == 0 && j >= last.from ) j--;
	
	/* synchronize alignments; note: there might be a gap at first.from or last.from, that's why we explicitely stop here  */
	while (ali->ali[i] != ali->ali[j] && i >= first.from && j >= last.from) {
	    while ( (ali->ali[i] > ali->ali[j]) ) { i --; s += ali->scores[i];}
	    while ( (ali->ali[j] > ali->ali[i]) ) { j --; s -= ali->scores[i];}
	}
	if ( i < first.from) i = first.from;
	if ( j < last.from)  j = last.from;

	if (verbose > LL2) 
	    printf ("%5.2f\t%5.2f\t%i\t%5i\t%5i\t%5i ", ali->scores[i], ali->scores[j], i, j, ali->ali[i], ali->ali[j]);
	
	if (ali->ali[i] == ali->ali[j])                 /* for this stupid border case such a fuss*/
	    s += ali->scores[i] - ali->scores[j];	/* calculate top - bottom */

	gapcost = gapvote[ali->ali[i]];
	s += lastgapcost - gapcost;			/* since gapcost is negative, "add" new gapcost and substract last gap-cost */
	if (verbose > LL2)
	    printf ("%5.2f\t%5.2f\t%5.2f\n", s, gapcost, lastgapcost);
	lastgapcost = gapcost;

	if (s > bestscore) {
	    bestscore = s;
	    bestcol  = i;
	}
	i--;
	j--;

    }

    free (gapvote);
    bestcol = ali->ali[bestcol];
    
    if (verbose > LL2) 
	printf("Chosen: bestcolumn %i with score %5.2f\n", bestcol, bestscore);
    
    return (bestcol);
} 

/*--------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
  Title: VoteGap
  explore region after $first_from and before $last_from and cut according to
  best gap
  The assumption is, that the repeats are more or less complete
  ---------------------------------------------------------------------------*/

int VoteGap ( first, last, ali, lrepeat ) 
SEQUENCEWINDOW first, last;
ARRAYALI *ali;
int lrepeat;
{
    CALCTYPE bestscore = 0;
    int bestcol   = first.from;
    CALCTYPE voted;
    int j = 0;
    int i = first.from;
    int gapregion;

    CALCTYPE * gapvote = GetGapVote( first.from, last.to, ali, lrepeat );

    if (verbose > LL2) 
	printf("Votegap called\ni\tj\tali_i\tgapvote\n");
  
    gapregion = (lrepeat * exploreratio);

    while (j < gapregion) {			/*explore region after start of first repeat */
	
	while ((ali->ali[i] == 0) && (i <= ali->length)) i++;

	voted = -gapvote[ali->ali[i]];          /* the more negative, the better, since the more gaps get annihilated */

	if (voted > bestscore) {
	    bestscore = voted;
	    bestcol   = i;
	}
	
	if (verbose > LL2) 
	    printf ("%i\t%i\t%i\t%5.2f\n", i, j, ali->ali[i], voted);

	i++;
	j++;
    }
    
    j = 0;
    i = last.from;
  
    while (j < gapregion) {			/*explore region before start of last repeat */
	while ((ali->ali[i] == 0) && (i > 0)) {
	    i--;
	}
	voted = -gapvote[ali->ali[i]];
	if (voted > bestscore) {
	    bestscore = voted;
	    bestcol = i;
	}
	
	if (verbose > LL2) 
	    printf ("%i\t%i\t%i\t%5.2f\n", i, j, ali->ali[i], voted);

	i--;
	j++;
    }
    
    free ( gapvote );
    
    if (verbose > LL2) 
	printf("Selected: bestcol %i with score %5.2f\n", bestcol, bestscore);

    return (bestcol);
}

/*--------------------------------------------------------------------
  explores the regions at the ends of the repeat-region, assuming that both
  ends are noisy and therefore tries to maximize the score of the internal
  repeat structure by pushing in/pulling out residues at the end.
  Since the alignment of internal residues never changes, the overall score
  differs just by the difference in the score of the residue pushed in and 
  the one popped out.
  Start: beginning of the first repeat is not part of the alignment is pulled
  in.
*/

int VotePushGap (first, last, ali, lrepeat) 
SEQUENCEWINDOW first;
SEQUENCEWINDOW last;
ARRAYALI *ali;
int lrepeat;
{

    /* score maximum columns against profile # note: approximation */
    /* move forward the boundary, and remember maximum score */
    int i = first.to;
    int j = last.to;
    CALCTYPE s = 0;
    int bestcol          = first.to + 1;	/* start: non-overlapping alignment; careful, might be gap */
    CALCTYPE bestscore   = 0;
    CALCTYPE lastgapcost = 0;           
    CALCTYPE gapcost;

    CALCTYPE *gapvote = (CALCTYPE*)GetGapVote( first.from, last.to, ali, lrepeat);
    if (verbose > LL2) 
	printf("Case VPG-1:\n scorei\tscorej\ti\tj\tali_i\tali_j\ts\tgapcost\tlastgapcost\n");

    while (i >= first.from) { 

	/* calculate difference scores for residues */
	/* move through gaps */
	while (ali->ali[i] == 0 && i >= first.from ) i--;
	while (ali->ali[j] == 0 && j >= last.from ) j--;
	
	/* synchronize alignments; note: there might be a gap at first.from or last.from, that's why we explicitely stop here  */
	while (ali->ali[i] != ali->ali[j] && i >= first.from && j >= last.from) {
	    while ( (ali->ali[i] > ali->ali[j]) ) { i --; s += ali->scores[i];}
	    while ( (ali->ali[j] > ali->ali[i]) ) { j --; s -= ali->scores[i];}
	}
	if ( i < first.from) i = first.from;
	if ( j < last.from)  j = last.from;

	if (verbose > LL2) 
	    printf ("%5.2f\t%5.2f\t%i\t%5i\t%5i\t%5i ", ali->scores[i], ali->scores[j], i, j, ali->ali[i], ali->ali[j]);
	
	if (ali->ali[i] == ali->ali[j])                 /* for this stupid border case such a fuss*/
	    s += ali->scores[i] - ali->scores[j];	/* calculate top - bottom */

	gapcost = gapvote[ali->ali[i]];
	s += lastgapcost - gapcost;			/* since gapcost is negative, "add" new gapcost and substract last gap-cost */
	if (verbose > LL2)
	    printf ("%5.2f\t%5.2f\t%5.2f\n", s, gapcost, lastgapcost);
	lastgapcost = gapcost;

	if (s > bestscore) {
	    bestscore = s;
	    bestcol  = i;
	}
	i--;
	j--;

    }

    free (gapvote);
    
    if (verbose > LL2) 
	printf("Chosen: bestcolumn %i with score %5.2f\n", bestcol, bestscore);
    
    return (bestcol);
} 


/*--------------------------------------------------------------------------*/
/*  Assume, that the alignment starts at the first residue of the first */
/*  repeat. Now start pulling the first repeat out and score each new */
/*  state. Choose the one with the minimum penalty/maximum score */
/*  The score is calculated as follows: */

/*  initially the score is the - gap-penalty of first residue */

/*  pull out one residue of first repeat   -> Penalty: The profile-score is subtracted. */
/*  add one residue of last repeat         -> Bonus:   The profile-score is added. */
/*  consider gap-preference of new start   -> Bonus:  subtract gap-penalty for new start */
/*  reconsider gap-preference of old start -> Pentalty: add gap-penalty */

/*  continue, until no residue is left from the last repeat. */
/*  to make things easier to write and more difficult to understand, */
/*  I will do it in reverse. */
/*-----------------------------------------------------------------------------*/

CALCTYPE PushGapVote ( ali, mincol, maxcol, lsequence, repeatscores, nrepeats ) 
ARRAYALI *ali;
int mincol, maxcol, lsequence;
CALCTYPE *repeatscores;
int *nrepeats;
{
    
    int lrepeat = maxcol - mincol + 1;  /*  length of repeat = length of profile */
    ARRAYALI *tempali, *shiftedali;
    int startcol_vpg, startcol_gap;
    CALCTYPE bestscore_bw, bestscore_vpg, bestscore_gap;
    int i;
    
    
    SEQUENCEWINDOW first, last, ol_first;
    verbose = LL1;

    startcol_vpg = startcol_gap = mincol;

    /* let alignment start at beginning of first repeat */
    first   = GetFirstRepeat( ali );
    shiftedali = (ARRAYALI*)CopyArrayAli( ali ); 

    ShiftAli ( ali->ali[first.from], 
	       mincol, maxcol, 
	       shiftedali, 
	       repeatscores,
	       nrepeats);
    
    if (verbose > LL3) {
	DumpArrayAlignment( shiftedali, "Temporary alignment" );
	printf ("Scores of %i repeats: \n", *nrepeats + 1);
	for (i = 0; i < *nrepeats; i++) printf("%5.2f ", repeatscores[i]);
	printf("\n");
    }

    if (*nrepeats <= 2) {
	if (verbose > LL1) 
	    printf ("Just two repeats, no register calculations necessary\n");
	FreeArrayAliMemory( shiftedali );
	bestscore_bw = ShiftAli( mincol, 
				 mincol, maxcol, 
				 ali, 
				 repeatscores, 
				 nrepeats);
	return (bestscore_bw);
    }
    

    /* get boundaries of first and last repeat */
    first  = GetFirstRepeat( shiftedali );    
    last   = GetLastRepeat ( shiftedali ); 
    if (verbose > LL1) 
	printf ("Range of first repeat %i-%i\nRange of last repeat: %i-%i\n", first.from, first.to, last.from, last.to);
	
    /* reduce to overlap region */
    ol_first.to = first.to;
    while ( (shiftedali->ali[last.to] < shiftedali->ali[ol_first.to]) || 
	    (shiftedali->ali[ol_first.to] < 1) 
	    ) {
	ol_first.to --;
    }
    
    if (verbose > LL1)
	printf ("After reduction: %i-%i and %i-%i\n", first.from, ol_first.to, last.from, last.to);
    
    /* do all shifts and take the best, this is almost exhaustive */

    startcol_gap = VoteGap( first, last, shiftedali, lrepeat);   
    if (startcol_gap != 0) {
	startcol_gap = ali->ali[startcol_gap];
	tempali   = (ARRAYALI*)CopyArrayAli( ali ); /* another working copy */
	bestscore_gap = ShiftAli( startcol_gap, 
				  mincol, maxcol, 
				  tempali, 
				  repeatscores, 
				  nrepeats);
	FreeArrayAliMemory( tempali );
    }
    
    first.to = ol_first.to; 
    startcol_vpg = VotePushGap( first, last, shiftedali, lrepeat); 
    startcol_vpg = ali->ali[startcol_vpg];
    if (startcol_vpg != 0) {
	tempali   = (ARRAYALI*)CopyArrayAli( ali ); /* another working copy */
	bestscore_vpg = ShiftAli( startcol_vpg, 
				  mincol, maxcol, 
				  tempali, 
				  repeatscores, 
				  nrepeats);
	FreeArrayAliMemory( tempali );
    }

    tempali = (ARRAYALI*)CopyArrayAli( ali ); /* another working copy */
    bestscore_bw  = ShiftAli( mincol, 
			      mincol, maxcol, 
			      tempali, 
			      repeatscores, 
			      nrepeats);
    FreeArrayAliMemory( tempali );
    
    if (verbose > LL0) 
	printf("Register: comparing scores: bestwindow = %5.2f (%i); gap = %5.2f (%i); vpg = %5.2f (%i)\n", 
	       bestscore_bw, mincol, bestscore_gap, startcol_gap, bestscore_vpg, startcol_vpg);
    
    if (bestscore_bw > bestscore_gap)  { startcol_gap  = mincol;       bestscore_gap = bestscore_bw; }
    if (bestscore_vpg > bestscore_gap) { startcol_gap  = startcol_vpg; bestscore_gap = bestscore_vpg;}
    
    FreeArrayAliMemory( shiftedali );

    /* now do the final shift */
    ShiftAli( startcol_gap, 
	      mincol, maxcol, 
	      ali, 
	      repeatscores, 
	      nrepeats);

    return (bestscore_gap);
}

/* 
   insertions in profile: from the point of view of bestwindow
   ---...------
   ---...------
   ---...------
   ------------ <- bestwindow
      +++

   insertions in sequence: from the point of view of bestwindow
   -----....-------
   -----....-------
   ----------------
   -----....------- <- bestwindow
            +

the columns marked + get a postive vote
j: gap penalties in sequence
i: gap penalties in profile

Note: since gap penalties are negative, the more negative the vote, the better.

Note: again there are two different objective functions to fulfill, depending on
whether the resulting overhanging ends are regarded as incomplete repeats or 
noise.

A. They are incomplete repeats. Then cuts are only allowed in the region, where
the first and last repeat not overlap.
*/

CALCTYPE UltimateRegister ( ali, mincol, maxcol, lsequence, repeatscores, nrepeats ) 
ARRAYALI *ali;
int mincol, maxcol, lsequence;
CALCTYPE *repeatscores;
int *nrepeats;
{
    CALCTYPE bestscore;
    int bestcol, newcol, lastcol;
    int i,j;
    int ngaps;
    int offset   = mincol - 1; 		/* for remapping from 1 to lprofile */
    int lprofile = maxcol - mincol + 1;
    CALCTYPE * profilepenalties = Malloc( (lprofile + 1) * sizeof(CALCTYPE ));

    for (i = 0; i <= lprofile; i++) profilepenalties[i] = 0;

    /* advance to first aligned residue */
    i = 1; while ( i <= ali->length && ali->ali[i] == 0) i++;

    lastcol = lprofile;
    newcol = ali->ali[i] - offset;
    while ( i <= ali->length) {

	/* calculate for insertions in profile */
	if (newcol < lastcol) { 	/* wrapping around case, insertions in profile */
	    if ( lprofile - lastcol > 0) 
		for (j = lastcol + 1; j <= lprofile; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( lprofile - lastcol) ;
	    if (newcol > 1)
		for (j = 1; j < newcol; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( newcol - 1) ;
	} else {                           
	    if (newcol > lastcol + 1 ) { /* insertion in profile */
		for (j = lastcol + 1; j < newcol; j++)
		    profilepenalties[j] += rg_gapopen_i + rg_gapelon_i * ( newcol - lastcol - 1 ) ;
	    }
	}
	
	lastcol = newcol;
	i++;
	ngaps = 0; while ( i <= ali->length && ali->ali[i] == 0) { i++; ngaps ++; }
	if (i > ali->length)
	    break;

	newcol = ali->ali[i] - offset;

	if (ngaps != 0) {
	    profilepenalties[newcol] += rg_gapopen_j + rg_gapelon_j * ( ngaps );
	}
    }

    if (verbose > LL1) 
	for (i = 1; i <= lprofile; i++) printf("%i (%i) %5.2f\n", i, i + offset, profilepenalties[i]);

    /* find optimum cutting point, i.e. the lowest score */
    bestcol = 1;
    bestscore = 0;
    for (i = 1; i <= lprofile; i++) 
	if (bestscore > profilepenalties[i]) {
	    bestscore = profilepenalties[i];
	    bestcol   = i;
	}

    if (verbose > LL0) 
	printf("Selected %i (%i) with score %5.2f\n", bestcol, bestcol + offset, bestscore);
    free( profilepenalties );
    
    
    
    bestscore = ShiftAli( bestcol + offset, 
			  mincol, maxcol, 
			  ali, 
			  repeatscores, 
			  nrepeats);



    return (bestscore);

}

CALCTYPE NewRegister ( ali, mincol, maxcol, lsequence, repeatscores, nrepeats ) 
ARRAYALI *ali;
int mincol, maxcol, lsequence;
CALCTYPE *repeatscores;
int *nrepeats;
{
    
    int lrepeat = maxcol - mincol + 1;  /*  length of repeat = length of profile */
    CALCTYPE bestscore;
    int bestcol;
    ARRAYALI *shiftedali;
    int startcol;
    int i;
    SEQUENCEWINDOW first, last, ol_first, region;

    startcol = mincol;
    /* let alignment start at beginning of first repeat */
    first   = GetFirstRepeat( ali );
    shiftedali = (ARRAYALI*)CopyArrayAli( ali ); 
    
    bestscore = ShiftAli ( ali->ali[first.from], 
			   mincol, maxcol, 
			   shiftedali, 
			   repeatscores,
			   nrepeats);
    
    if (verbose > LL2) {
	DumpArrayAlignment( shiftedali, "Temporary alignment" );
	printf ("Scores of %i repeats: \n", *nrepeats + 1);
	for (i = 0; i < *nrepeats; i++) printf("%5.2f ", repeatscores[i]);
	printf("\n");
    }
    
    if (*nrepeats <= 2) {
	if (verbose > LL1) 
	    printf ("Just two repeats after shifting, register starts at first repeat\n");
	FreeArrayAliMemory( shiftedali );
	bestscore = ShiftAli( ali->ali[first.from], 
			      mincol, maxcol, 
			      ali, 
			      repeatscores, 
			      nrepeats);
	return (bestscore);
    }
    
    startcol = mincol;

    /* get boundaries of first and last repeat */
    first  = GetFirstRepeat( shiftedali );    
    last   = GetLastRepeat ( shiftedali ); 
    if (verbose > LL2) 
	printf ("Range of first repeat %i-%i\nRange of last repeat: %i-%i\n", first.from, first.to, last.from, last.to);
	
    /* reduce to overlap region */
    ol_first.to = first.to;
    while ( (shiftedali->ali[last.to] < shiftedali->ali[ol_first.to]) || 
	    (shiftedali->ali[ol_first.to] < 1) 
	    ) {
	ol_first.to --;
    }
    
    if (verbose > LL2)
	printf ("After reduction: %i-%i and %i-%i\n", first.from, ol_first.to, last.from, last.to);
    
    if (                                            /*  if the repeats are mostly complete */
	(shiftedali->ali[last.to] * compfactor >= lrepeat )
	) {              
	if (verbose > LL1) printf("Cutting in region D (non-overlap): overlap-length %i lrepeat %i\n", shiftedali->ali[last.to], lrepeat);
	   
	region.from = shiftedali->ali[ol_first.to]; /* ends of non-overlapping region (well, slight overlap of one residue) */
	region.to   = lrepeat;
	bestcol     = CutInD( region, shiftedali, lrepeat );	/* push ends into alignment */
	
	/* note: if there are masked repeats, the cut might have been placed at a profile-position, that does not exist, therefore remember next smallest position */
/* 	i = mincol;  */
	/* 	lastcol  = mincol; */
/* 	startcol++; */
/* 	while (i <= maxcol) { */
/* 		if (thiscol = shiftedali->ali[startcol] == bestcol) { */
/* 		    startcol = i; */
/* 		    break; */
/* 		if (thiscol < bestcol && thiscol > lastcol) { */
/* 		    lastcol  = thiscol; */
/* 		    startcol = i; */
/* 		} */
/* 	    } */
 	     
    } else {                                        
	if (verbose > LL1) printf("Cutting in region E (overlap): overlap-length %i lrepeat %i\n", shiftedali->ali[last.to], lrepeat);
	first.to = ol_first.to;                     /* set first to end of overlapping region */
	bestcol  = CutInE( first, last, shiftedali, lrepeat);	/*  deal with noisy ends */
    }

    startcol = mincol; while (startcol <= maxcol && shiftedali->ali[startcol] != bestcol) startcol++; /* find starting column in alignment */
    
    /* but: there are cases, when the best window is not found, in this case set startcol to mincol */
    if (startcol == 0) {
	startcol = mincol;
	if (verbose > LL1) 
	    printf ("No residue in bestwindow, setting startcol to mincol\n");
    }
    
    if (verbose > LL2)
        printf ("Startcolumn for shifting alignment: %i\n", startcol );
    
    /* now do the final shift */

    bestscore = ShiftAli( startcol,
                          mincol, maxcol,
                          ali,
                          repeatscores,
                          nrepeats);
    
    
    FreeArrayAliMemory( shiftedali ); 
    return (bestscore);  

}

/*--------------------------------------------------------------------------*/
/* Title  : CountRepeats                                                   */
/* Funct. : count repeats in a profile-alignment                         */
/* Author : Andreas Heger		                               */
/* Created: 9.1.1999                                                    */
/*--------------------------------------------------------------------- */	    

int CountRepeats( ali ) 
    ARRAYALI *ali;
{
    
    int i, lastres;
    int nrepeats = 0;

    int lali = ali->length;
    
    /* go to beginning of first repeat */
    i = 1; while ( (i <= lali) && (ali->ali[i] == 0) ) i++;
    if (i > lali) return (0);

    /* intialize */
    lastres  = 0;
    
    /* start loop */
    while (i <= lali) {
	
	while ((i <= lali) && ali->ali[i] == 0) i++; 	           /* count gaps and advance to next non-gap character */
	if (i > lali) break;
	
	if (ali->ali[i] < lastres ||				/* end previous repeat and start a new repeat */
	    ( (lastres - ali->ali[i]) > splitdistance) )            
	    nrepeats++;
	
	lastres  = ali->ali[i];                                        
	i++;                                                        /* advance to next residue */
    }
    
    return (nrepeats+1);
}
/*--------------------------------------------------------------------------*/
/* Title  : CountAndScoreRepeats                                        */
/* Funct. : count and scores the repeats in a profile-alignment         */
/* Author : Andreas Heger		                               */
/* Created: 9.1.1999                                                    */
/*--------------------------------------------------------------------- */	    

int CountAndScoreRepeats( ali, lrepeat, gapopen, gapelon, windows ) 
    ARRAYALI *ali;
int lrepeat;
CALCTYPE gapopen, gapelon;
SEQUENCEWINDOW *windows;
{

    int i;
    CALCTYPE s;
    int lastres, newres, ngaps, lasti;
    int nrepeats = 0;
    int lali = ali->length;
    
    /* go to beginning of first repeat */
    i = 1; while ( (i <= lali) && (ali->ali[i] == 0) ) i++;
    if (i > lali) return (0);

    /* intialize */
    lastres  = 0;
    newres = ali->ali[i];

    /* compute gap penalty for first repeat */
    if (newres > 1) 
	s = gapopen + gapelon * (newres - 1);
    else 
	s = 0;

    windows[nrepeats].from = i;
    
    /* start main loop */
    while (i <= lali) {
	
	ngaps = 0;
	while ((i <= lali) && ali->ali[i] == 0) { 	           /* count gaps and advance to next non-gap character */
	    ngaps++;
	    i++;
	}
	if (i > lali) break;
	
	newres = ali->ali[i];
	
	if (newres < lastres || ngaps > splitdistance) {	   /* end previous repeat and start a new repeat, if wrapped over register or big gap */

	    if (lrepeat - lastres > 0) {	                   /* penalize for deletions in profile at the end of the repeat, only if wrapped */
		s += gapopen + gapelon * (lrepeat - lastres );
	    }
	    /* printf(" Added new repeat %i %i %i %i %i %i\n", i, lasti, ngaps, newres, lastres, nrepeats); */
	    windows[nrepeats].to   = lasti;
	    windows[nrepeats].score= s;                             /* save score */
	    nrepeats++;
	    windows[nrepeats].from = i;

	    if (newres > 1) 	                                    /* compute gap penalty for start of repeat */
		s = gapopen + gapelon * (newres - 1);
	    else 
		s = 0;
	} else {                                                    /* else penalize for gaps if continuation of same repeat*/
	    if (ngaps) s += gapopen + ngaps * gapelon;
	}

	s += ali->scores[i];                                        /* add alignment score to repeat score */
	lasti    = i;
	lastres  = newres;                                        
	i++;                                                        /* advance to next residue */
    }

    if (lrepeat - lastres > 0) 
	s += gapopen + gapelon * (lrepeat - lastres );

    windows[nrepeats].to = lasti;
    windows[nrepeats].score = s;

    if (verbose > LL2) { 
	printf("Repeats with score from-to; score\n");
	for (i = 0; i <= nrepeats; i++) 
	    printf("%i-%i; %5.2f\n", windows[i].from, windows[i].to, windows[i].score);
	printf("\n");
    }

    return (nrepeats+1);
}


/*--------------------------------------------------------------------------*/
/* Title  : CountAndScoreRepeats                                        */
/* Funct. : Shifts a profile-alignment, so that the column startcol is  */
/*          mapped to 1                                                 */
/* Author : Andreas Heger		                               */
/* Created: 9.1.1999                                                    */
/*--------------------------------------------------------------------- */	    
void ShiftProfileAli ( startcol, lprofile, ali )
int startcol, lprofile;
ARRAYALI *ali;
{
    int i, j;
    
    int * map = (int*)Malloc((lprofile + 1) * sizeof(int));

    /* build map */
    j = 0;
    for (i = 0; i <= lprofile; i++) map[i] = 0;
    for (i = startcol; i <= lprofile; i++) map[i]=++j; 
    for (i = 1; i <  startcol; i++)        map[i]=++j; 
    
    /* do the mapping */
    i = 1;
    while (i <= ali->length) {
	ali->ali[i] = map[ali->ali[i]];
	i++;
    }
    free (map);
}    

/* get a profile-alignment as input */
void NewestRegister ( ali, lrepeat, anchorres ) 
ARRAYALI *ali;
int lrepeat;
int anchorres;
{
    int startcol;
    SEQUENCEWINDOW first, last, ol_first, region;
    int nrepeats;

    /* let alignment start at beginning of first repeat */
    first   = GetFirstRepeat( ali );
    ShiftProfileAli ( ali->ali[first.from], lrepeat, ali );
    
    if (verbose > LL2) 
	DumpArrayAlignment( ali, "Alignment on first repeat" );
    
    nrepeats = CountRepeats( ali );
    
    if (nrepeats <= 2) {
	if (verbose > LL1) 
	    printf ("Just two repeats after shifting, register starts at first repeat\n");
	return;
    }
    
    startcol = 1;

    /* get boundaries of first and last repeat */
    first  = GetFirstRepeat( ali );    
    last   = GetLastRepeat ( ali ); 
    if (verbose > LL2) 
	printf ("Range of first repeat %i-%i\nRange of last repeat: %i-%i\n", first.from, first.to, last.from, last.to);
	
    /* reduce to overlap region */
    ol_first.to = first.to;
    while ( (ali->ali[last.to] < ali->ali[ol_first.to]) || 
	    (ali->ali[ol_first.to] < 1) 
	    ) {
	ol_first.to --;
    }
    
    if (verbose > LL2)
	printf ("After reduction: %i-%i and %i-%i\n", first.from, ol_first.to, last.from, last.to);
    
    if (             
	(ali->ali[ol_first.to] == lrepeat)) {
	/* no register calculation, if repeats are complete (is not handeled well by CutInD) */
	if (verbose > LL1) printf("Complete overlap, no register calculation necessary\n");
	startcol = 1;
    } else if (ali->ali[last.to] * compfactor >= lrepeat ) {
	/* if the repeats are mostly complete cut in non-overlapping region */
	if (verbose > LL1) printf("Cutting in region D (non-overlap): overlap-length %i lrepeat %i\n", ali->ali[last.to], lrepeat);

	region.from  = ali->ali[last.to];			/* ends of non-overlapping region */
	region.to    = lrepeat;
	startcol     = CutInD( region, ali, lrepeat );		/* push ends into alignment */
	if (startcol == 0) startcol = 1;                        /* choose start of first repeat, if no decision */
    } else {                                        
	/* if the repeats are overlap only slightly, cut in overlapping region */
	if (verbose > LL1) printf("Cutting in region E (overlap): overlap-length %i lrepeat %i\n",     ali->ali[last.to], lrepeat);
	first.to = ol_first.to;                                 /* set first to end of overlapping region */
	startcol  = CutInE( first, last, ali, lrepeat);		/*  deal with noisy ends */
	if (startcol == 0) startcol = ali->ali[anchorres];               /* start at bestwindow, if no decision */
    }
    
    /* but: there are cases, when the bestwindow is not found, in this case do not change anything */
    if (startcol == 0) {
	if (verbose > LL1) 
	    printf ("No residue in bestwindow, no change\n");
	startcol = ali->ali[anchorres];
    }
    
    if (verbose > LL2)
        printf ("Startcolumn for shifting alignment: %i\n", startcol );
    
    /* now do the final shift */
    
    ShiftProfileAli( startcol, lrepeat,	ali);
    
}
