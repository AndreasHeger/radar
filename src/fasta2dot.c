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
$Id: fasta2dot.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/
#include "align.h"
#include "toolbox.h"
#include "align_tools.h"
#include <stdio.h>
#include <stdlib.h>


#define MAX(x,y) (x > y) ? x : y
#define MYINDEX(x,y) x * M + y

#define MYINDEXSAVESPACE(x,y) x * M - (int)((x-1) * x / 2) + y - x
/* #define MYINDEXSAVESPACE(x,y) x  * ( M - ( x - 1 ) / 2 ) + y : does not work, because of type-conversion*/

int GetStart( buffer, ali ) 
    char * buffer ;
char * ali;

{
    int i = 3;
    int imax = strlen( buffer );
    int number = 0;
    int lastoccurance = 0;
    int ngaps = 0;

    if (imax < 5) /* if there are no line numbers in this line */
	return 0;

    while( i < imax) {
	if (buffer[i] >= 48 && buffer[i] <= 57 ) {
	    number = number * 10 + (buffer[i] - 48);
	    lastoccurance = i;
	} else {
	    if (lastoccurance) break;
	}
	i++;
    }

    if (i == imax)       /* no number in this line (happens with short fragments */ 
	return 0;

    if (number == 0) 
	return 0;
    /* correct for gaps in the first part of the sequence*/
    for (i = 7; i < lastoccurance; i++) {
	if (ali[i] == '-') ngaps++;
    }

    return (number - lastoccurance + 7 + ngaps);
}
 
int Fasta2Dots( row, col, score, source, profile, M, mindiag )
int *row, *col;
CALCTYPE *score;
char *source;
PROFILECOLUMN *profile;
int M;
int mindiag;

{
    FILE *fh;
    char buffer[MAXLENBUFFER];
    char buffer1[MAXLENBUFFER];
    char buffer2[MAXLENBUFFER];
    SEQTYPE *decode = GetDecode();
    CALCTYPE maxscore;
    
    char * matrix = (char*)Malloc (sizeof(char) * (M + 1)  * (M + 1)); /* quick hack, not space efficient */
    
    int r1, r2, start1, start2, i;
    int res1, res2;
    char s1, s2;
    int ndots = 0;

#ifdef DEBUG
    printf("Fasta2Dot called with mindiag %i:\n", mindiag);
    PrintProfile( profile, M);
#endif

    res1 = res2 = 0;
    
    for (i = 0; i < (M+1)*(M+1); i ++) matrix[i] = 0;

    OpenFile( &fh, source, "r", "fasta2dot");
    
    while (!feof(fh)) {
	fgets( buffer, MAXLENBUFFER, fh);
	/* found beginning of alignment */
	if (strstr( buffer, "identity")) {
	    fgets(buffer, MAXLENBUFFER, fh);             /* blank line */ 
	    fgets(buffer, MAXLENBUFFER, fh);             /* read top scale */

	    res1   = res2   = 0;                         /* reset counters, important!! */
	    start1 = start2 = 0;

	    while ( buffer[0] != '-') {
		fgets(buffer1, MAXLENBUFFER, fh); /* read first line of alignment */
		start1 = GetStart( buffer, buffer1 );
		if (start1 != 0) res1 = start1;   /* if just a short fragment, carry over from last alignment */
		fgets(buffer2, MAXLENBUFFER, fh); 
		fgets(buffer2, MAXLENBUFFER, fh); /* read second line of alignment */
		fgets(buffer,  MAXLENBUFFER, fh); /* read bottom scale */
		start2 = GetStart( buffer, buffer2 );
		if (start2 != 0) res2 = start2;   /* if just a short fragment, carry over from last alignment */
		
#ifdef DEBUG
		printf("Parsing: start1 %i start2 %i\n", res1, res2);
#endif
		/* safety check for very, very short alignments, where startresidues can not be found, can occur */
		if (res1 < 1 || res2 < 1 || res1 > M || res2 > M) 
		    break; /* skip */
		for (i = 7; i < strlen(buffer1); i++) {
		    s1 = buffer1[i];
		    if ((s1 < 65 || s1 > 90) && s1 != 45) break;
		    s2 = buffer2[i];
		    if ( (abs(res2 - res1) > mindiag) &&
			 s1 != '-' && 
			 s2 != '-') {
			r1 = decode[s1]; r2 = decode[s2];
			maxscore = MAX(profile[res1][r2], profile[res2][r1]); 
#ifdef DEBUG
			printf ("-%i-%i-%c-%c-%i-%i-%5.2f %5.2f %5.2f\n", res1, res2, s1, s2, r1, r2, maxscore, profile[res1][r2], profile[res2][r1]); 
#endif

			if (res1 < res2) { /* only use upper diagonal in matrix*/
			    r1 = res1; 
			    r2 = res2;
			} else {
			    r1 = res2;
			    r2 = res1;
			}
			
			if (matrix[MYINDEX(r1, r2)] == 0) {
			    row[ndots]   = res1; 
			    col[ndots]   = res2;
			    score[ndots] = maxscore;
			    ndots++;
			    
			    row[ndots]   = res2; 
			    col[ndots]   = res1;
			    score[ndots] = maxscore;
			    ndots ++;
			    matrix[MYINDEX(r1,r2)] = 1;
			}
		    }
		    
		    if (buffer1[i] != '-') res1++;
		    if (buffer2[i] != '-') res2++;
		}
		fgets(buffer,  MAXLENBUFFER, fh); /* blank line */
		fgets(buffer,  MAXLENBUFFER, fh); /* read top scale or ------ */
	    }
	}
    }

    FClose( fh );
    free (matrix );

#ifdef DEBUG
    PrintDots( row, col, score, ndots);
#endif
    

    return ndots;
}

/* score Profile-Columns against each other */
int Fasta2DotsCorrect( row, col, score, source, profile, counts, M, mindiag )
int *row, *col;
CALCTYPE *score;
char *source;
PROFILECOLUMN *profile;
COUNTCOLUMN *counts;
int M;
int mindiag;

{
    FILE *fh;
    char buffer[MAXLENBUFFER];
    char buffer1[MAXLENBUFFER];
    char buffer2[MAXLENBUFFER];
    SEQTYPE *decode = GetDecode();
    CALCTYPE maxscore;
    FREQUENCYCOLUMN * frequencies = (FREQUENCYCOLUMN*)Malloc( sizeof(FREQUENCYCOLUMN) * (M + 1) );
    char * matrix = (char*)Malloc (sizeof(char) * (M + 1)  * (M + 1)); /* quick hack, not space efficient */

    int r1, r2, start1, start2, i;
    int res1, res2;
    char s1, s2;
    int ndots = 0;

#ifdef DEBUG
    printf("Fasta2DotsCorrect called with mindiag %i:\n", mindiag);
    PrintProfile( profile, M);
    PrintCounts( counts, M );
#endif

    Counts2Frequencies( counts, M, frequencies);
#ifdef DEBUG
    PrintFrequencies( frequencies, M);
#endif
    
    for (i = 0; i < (M+1)*(M+1); i ++) matrix[i] = 0;

    OpenFile( &fh, source, "r", "fasta2dot");
    
    while (!feof(fh)) {
	fgets( buffer, MAXLENBUFFER, fh);
	/* found beginning of alignment */
	if (strstr( buffer, "identity")) {
	    fgets(buffer, MAXLENBUFFER, fh);             /* blank line */ 
	    fgets(buffer, MAXLENBUFFER, fh);             /* read top scale */
	    
	    res1   = res2   = 0;                         /* reset counters, important!! */
	    start1 = start2 = 0;

	    while ( buffer[0] != '-') {
		fgets(buffer1, MAXLENBUFFER, fh); /* read first line of alignment */
		start1 = GetStart( buffer, buffer1 );
		if (start1 != 0) res1 = start1;   /* if just a short fragment, carry over from last alignment */
		fgets(buffer2, MAXLENBUFFER, fh); 
		fgets(buffer2, MAXLENBUFFER, fh); /* read second line of alignment */
		fgets(buffer,  MAXLENBUFFER, fh); /* read bottom scale */
		start2 = GetStart( buffer, buffer2 );
		if (start2 != 0) res2 = start2;   /* if just a short fragment, carry over from last alignment */
		
#ifdef DEBUG
		printf("Parsing: start1 %i start2 %i\n", res1, res2);
#endif
		/* safety check for very, very short alignments, where startresidues can not be found, can occur */
		if (res1 < 1 || res2 < 1 || res1 > M || res2 > M) 
		    break; /* skip */
		for (i = 7; i < strlen(buffer1); i++) {
		    s1 = buffer1[i];
		    if ((s1 < 65 || s1 > 90) && s1 != 45) break;
		    s2 = buffer2[i];
		    if ( (abs(res2 - res1) > mindiag) &&
			 s1 != '-' && 
			 s2 != '-') {
			r1 = decode[s1]; r2 = decode[s2];
			maxscore = ProfileScore(profile[res1], frequencies[res1], profile[res2], frequencies[res2]); 
#ifdef DEBUG
			printf ("-%i-%i-%c-%c-%i-%i-%5.2f %5.2f %5.2f\n", res1, res2, s1, s2, r1, r2, maxscore, profile[res1][r2], profile[res2][r1]); 
#endif

			if (res1 < res2) { /* only use upper diagonal in matrix*/
			    r1 = res1; 
			    r2 = res2;
			} else {
			    r1 = res2;
			    r2 = res1;
			}

			if (matrix[MYINDEX(r1, r2)] == 0) {
			    row[ndots]   = res1; 
			    col[ndots]   = res2;
			    score[ndots] = maxscore;
			    ndots++;
			    
			    row[ndots]   = res2; 
			    col[ndots]   = res1;
			    score[ndots] = maxscore;
			    ndots ++;
			    matrix[MYINDEX(r1,r2)] = 1;
			}
		    }
		    
		    if (buffer1[i] != '-') res1++;
		    if (buffer2[i] != '-') res2++;
		}
		fgets(buffer,  MAXLENBUFFER, fh); /* blank line */
		fgets(buffer,  MAXLENBUFFER, fh); /* read top scale or ------ */
	    }
	}
    }

    FClose( fh );
    free (matrix );

#ifdef DEBUG
    PrintDots( row, col, score, ndots);
#endif
    
    free( frequencies );

    return ndots;
}


int Fasta2Dotfile( source, dest1, profile, M, mindiag )
    char *source;
char *dest1;
PROFILECOLUMN *profile;
int M;
int mindiag;

{
    int *row   = (int*)Malloc (MAXNDOTS * sizeof( int )); 
    int *col   = (int*)Malloc (MAXNDOTS * sizeof( int ));
    CALCTYPE *score = (CALCTYPE*)Malloc (MAXNDOTS * sizeof( CALCTYPE ));

    int ndots = Fasta2Dots( row, col, score, source, profile, M , mindiag);
    
    WriteDots( row, col, score, ndots, dest1 );

    free( row );
    free( col );
    free( score );

    return ndots;
}

/*-----------------------*/
void FillMatrix( matrix, M, source, mindiag ) 
    char * matrix;
int M;
char * source;
int mindiag;
{
    int r1, r2, start1, start2, i;
    int res1, res2;
    char s1, s2;
    FILE *fh;
    char buffer[MAXLENBUFFER];
    char buffer1[MAXLENBUFFER];
    char buffer2[MAXLENBUFFER];
    SEQTYPE *decode = GetDecode();

    OpenFile( &fh, source, "r", "fasta2dot");
    
    while (!feof(fh)) {
	fgets( buffer, MAXLENBUFFER, fh);
	/* found beginning of alignment */
	if (strstr( buffer, "identity")) {
	    fgets(buffer, MAXLENBUFFER, fh);             /* blank line */ 
	    fgets(buffer, MAXLENBUFFER, fh);             /* read top scale */
	    
	    res1   = res2   = 0;                         /* reset counters, important!! */
	    start1 = start2 = 0;

	    while ( buffer[0] != '-') {
		fgets(buffer1, MAXLENBUFFER, fh); /* read first line of alignment */
		start1 = GetStart( buffer, buffer1 );
		if (start1 != 0) res1 = start1;   /* if just a short fragment, carry over from last alignment */
		fgets(buffer2, MAXLENBUFFER, fh); 
		fgets(buffer2, MAXLENBUFFER, fh); /* read second line of alignment */
		fgets(buffer,  MAXLENBUFFER, fh); /* read bottom scale */
		start2 = GetStart( buffer, buffer2 );
		if (start2 != 0) res2 = start2;   /* if just a short fragment, carry over from last alignment */
#ifdef DEBUG
		printf("Parsing: start1 %i start2 %i\n", res1, res2);
#endif
		/* safety check for very, very short alignments, where startresidues can not be found, can occur */
		if (res1 < 1 || res2 < 1 || res1 > M || res2 > M) 
		    break; /* skip */
		for (i = 7; i < strlen(buffer1); i++) {
		    s1 = buffer1[i];
		    if ((s1 < 65 || s1 > 90) && s1 != 45) break;
		    s2 = buffer2[i];
		    if ( (abs(res2 - res1) > mindiag) &&
			 s1 != '-' && 
			 s2 != '-') {

			if (res1 < res2) { /* only use upper diagonal in matrix*/
			    r1 = res1; 
			    r2 = res2;
			} else {
			    r1 = res2;
			    r2 = res1;
			}

#ifdef DEBUG
			printf ("-%i-%i-%c-%c-%i-%i \n", res1, res2, r1, r2); 
#endif

			matrix[MYINDEX(r1,r2)] = 1;
		    }
		    
		    if (buffer1[i] != '-') res1++;
		    if (buffer2[i] != '-') res2++;
		}
		fgets(buffer,  MAXLENBUFFER, fh); /* blank line */
		fgets(buffer,  MAXLENBUFFER, fh); /* read top scale or ------ */
	    }
	}
    }
    FClose( fh );
}

/*------------------------------------------------------------------------------------------------*/
/* score Profile-Columns against each other */
int FastaFiles2DotsCorrect( row, col, score, source1, source2, profile, counts, M, mindiag )
int *row, *col;
CALCTYPE *score;
char *source1;
char *source2;
PROFILECOLUMN *profile;
COUNTCOLUMN *counts;
int M;
int mindiag;
{
    CALCTYPE maxscore;
    FREQUENCYCOLUMN * frequencies = (FREQUENCYCOLUMN*)Malloc( sizeof(FREQUENCYCOLUMN) * (M + 1) );

    char * matrix = (char*)Malloc (sizeof(char) * (M + 1)  * (M + 1)); /* quick hack, not space efficient */

    int ndots = 0;
    int i,j;


#ifdef DEBUG
    printf("FastaFiles2DotsCorrect called with mindiag %i:\n", mindiag);
    PrintProfile( profile, M);
    PrintCounts( counts, M );
#endif

    Counts2Frequencies( counts, M, frequencies);
#ifdef DEBUG
    PrintFrequencies( frequencies, M);
#endif
    
    for (i = 0; i < (M+1)*(M+1); i ++) matrix[i] = 0;

    FillMatrix( matrix, M, source1, mindiag) ;
    FillMatrix( matrix, M, source2, mindiag) ;

    for (i = 1; i < M; i++) {
	for ( j = i+1; j <= M; j++) {
	    if (matrix[MYINDEX(i, j)] == 1) {
		maxscore = ProfileScore(profile[i], frequencies[i], profile[j], frequencies[j]); 
		row[ndots]   = i; 
		col[ndots]   = j;
		score[ndots] = maxscore;
		ndots++;
			    
		row[ndots]   = j; 
		col[ndots]   = i;
		score[ndots] = maxscore;
		ndots ++;
	    }
	}
    }

    free (matrix );
    
#ifdef DEBUG
    PrintDots( row, col, score, ndots);
#endif
    
    free( frequencies );

    return ndots;
}

/*-----------------------*/
void FillMatrixSaveSpace( matrix, M, source, mindiag ) 
    char * matrix;
int M;
char * source;
int mindiag;
{
    int r1, r2, start1, start2, i;
    int res1, res2;
    char s1, s2;
    FILE *fh;
    char buffer[MAXLENBUFFER];
    char buffer1[MAXLENBUFFER];
    char buffer2[MAXLENBUFFER];
    SEQTYPE *decode = GetDecode();

    OpenFile( &fh, source, "r", "fasta2dot");
    
    while (!feof(fh)) {
	fgets( buffer, MAXLENBUFFER, fh);
	/* found beginning of alignment */
	if (strstr( buffer, "identity")) {
	    fgets(buffer, MAXLENBUFFER, fh);             /* blank line */ 
	    fgets(buffer, MAXLENBUFFER, fh);             /* read top scale */
	    
	    res1   = res2   = 0;                         /* reset counters, important!! */
	    start1 = start2 = 0;

	    while ( buffer[0] != '-' && !feof(fh)) {
		fgets(buffer1, MAXLENBUFFER, fh); /* read first line of alignment */
		start1 = GetStart( buffer, buffer1 );
		if (start1 != 0) res1 = start1;   /* if just a short fragment, carry over from last alignment */
		fgets(buffer2, MAXLENBUFFER, fh); 
		fgets(buffer2, MAXLENBUFFER, fh); /* read second line of alignment */
		fgets(buffer,  MAXLENBUFFER, fh); /* read bottom scale */
		start2 = GetStart( buffer, buffer2 );
		if (start2 != 0) res2 = start2;   /* if just a short fragment, carry over from last alignment */
		
#ifdef DEBUG 
		printf("Parsing: start1 %i start2 %i\n", res1, res2);
#endif 
		/* safety check for very, very short alignments, where startresidues can not be found, can occur */
		if (res1 < 1 || res2 < 1 || res1 > M || res2 > M) 
		    break; /* skip */
		for (i = 7; i < strlen(buffer1); i++) {
		    s1 = buffer1[i];
		    if ((s1 < 65 || s1 > 90) && s1 != 45) break;
		    s2 = buffer2[i];
		    if ( (abs(res2 - res1) > mindiag) &&
			 s1 != '-' && 
			 s2 != '-') {
			if (res1 < res2) { /* only use upper diagonal in matrix*/
			    r1 = res1 - 1; 
			    r2 = res2 - 1;
			} else {
			    r1 = res2 - 1;
			    r2 = res1 - 1;
			}
#ifdef DEBUG 
			printf ("-%i-%i-%c-%c-%i-%i\n", res1, res2, s1, s2, r1, r2); 
#endif
			/* printf ("%i %i %i\n", r1, r2, MYINDEXSAVESPACE(r1, r2));*/
			matrix[MYINDEXSAVESPACE(r1, r2)] = 1;
		    }
		    
		    if (buffer1[i] != '-') res1++;
		    if (buffer2[i] != '-') res2++;
		}
		fgets(buffer,  MAXLENBUFFER, fh); /* blank line */
		fgets(buffer,  MAXLENBUFFER, fh); /* read top scale or ------ */
	    }
	}
    }
    FClose( fh );
}


/*------------------------------------------------------------------------------------------------*/
/* score Profile-Columns against each other */
int FastaFiles2DotsCorrectSaveSpace( row, col, score, source1, source2, profile, counts, M, mindiag )
int *row, *col;
CALCTYPE *score;
char *source1;
char *source2;
PROFILECOLUMN *profile;
COUNTCOLUMN *counts;
int M;
int mindiag;
{
    CALCTYPE maxscore;
    FREQUENCYCOLUMN * frequencies = (FREQUENCYCOLUMN*)Malloc( sizeof(FREQUENCYCOLUMN) * (M + 1) );

    int matrixsize = M * (M + 1) / 2;
    char * matrix = (char*)Malloc (sizeof(char) * matrixsize); /* quick hack, not space efficient */

    int ndots = 0;
    int i,j;
    int i2, j2;

#ifdef DEBUG
    printf("FastaFiles2DotsCorrect called with mindiag %i:\n", mindiag);
    PrintProfile( profile, M);
    PrintCounts( counts, M );
    printf("Size of temporary matrix will be: %i\n", matrixsize);
#endif

    Counts2Frequencies( counts, M, frequencies);
#ifdef DEBUG
    PrintFrequencies( frequencies, M);
#endif
    
    for (i = 0; i < matrixsize; i ++) matrix[i] = 0;
    FillMatrixSaveSpace( matrix, M, source1, mindiag) ;
    FillMatrixSaveSpace( matrix, M, source2, mindiag) ;

    /* matrix: residue numbering is from 0 to M-1 */
    /* profile, etc.: residue numbering is from 1 to M */
     
    for (i = 0; i < M - 1; i++) {
	for ( j = i+1; j < M; j++) {
	    if (matrix[MYINDEXSAVESPACE(i, j)] == 1) {
		i2 = i+1;
		j2 = j+1;
		maxscore = ProfileScore(profile[i2], frequencies[i2], profile[j2], frequencies[j2]); 

		row[ndots]   = i2; 
		col[ndots]   = j2;
		score[ndots] = maxscore;
		ndots++;

		row[ndots]   = j2; 
		col[ndots]   = i2;
		score[ndots] = maxscore;
		ndots++;

#ifdef SAVE
		if (ndots == MAXNDOTS) {
		    printf( "Error: too much dots\n");
		    exit(0);
		}
#endif
/* 		row[ndots]   = j;  */
/* 		col[ndots]   = i; */
/* 		score[ndots] = maxscore; */
/* 		ndots ++; */
/* #ifdef SAVE */
/* 		if (ndots == MAXNDOTS) { */
/* 		    printf( "Error: too much dots\n"); */
/* 		    exit(0); */
/* 		} */
/* #endif */

	    }
	}
    }

    free (matrix );
    
#ifdef DEBUG
    PrintDots( row, col, score, ndots);
#endif
    
    free( frequencies );

    return ndots;
}

