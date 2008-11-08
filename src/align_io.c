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
$Id: align_io.c,v 1.3 2005/02/17 23:13:32 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#ifndef __ALIGN_H
#include "align.h"
#endif
 
#include "align_io.h"
#include "align_tools.h"
#include "toolbox.h"
#include "string.h"


/*--------------------------------------------------------------------*/
int ReadDots( filename, row, col, score ) 
char *filename;
int **row;
int **col;
CALCTYPE **score;

{
    FILE *fh;
    int M;
    long i;
    
    int *R;
    CALCTYPE *S;

    OpenFile( &fh, filename, "r", "align");
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of profile */
    M = (long)(i / ( 2 * sizeof(int) + sizeof(CALCTYPE)) );
    /* read in rows */ 
    R = (int *)Malloc( M * sizeof(int) ); 
    fread( R, sizeof( int ), M, fh);
    *row = R;
    /* read in columns */ 
    R = (int *)Malloc( M * sizeof(int) ); 
    fread( R, sizeof( int ), M, fh);
    *col = R;
    /* read in scores */ 
    S = (CALCTYPE *)Malloc( M * sizeof(CALCTYPE) ); 
    fread( S, sizeof( CALCTYPE ), M, fh);
    *score = S;

    fclose(fh);

    return M;
}
/*--------------------------------------------------------------------*/
int ReadDotsOnlyRowsCols( filename, row, col ) 
char *filename;
int **row;
int **col;

{
    FILE *fh;
    int M;
    long i;
    
    int *R;

    OpenFile( &fh, filename, "r", "align");
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of profile */
    M = (long)(i / ( 2 * sizeof(int) ));
    /* read in rows */ 
    R = (int *)Malloc( M * sizeof(int) ); 
    fread( R, sizeof( int ), M, fh);
    *row = R;
    /* read in columns */ 
    R = (int *)Malloc( M * sizeof(int) ); 
    fread( R, sizeof( int ), M, fh);
    *col = R;

    fclose(fh);

    return M;
}


/*--------------------------------------------------------------------*/
int ReadProfile( filename, A )
    char *filename;
    PROFILECOLUMN **A;
{
    FILE *fh;
    int M;
    long i;
    
    PROFILECOLUMN *R;

    OpenFile( &fh, filename, "r", "align");
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of profile */
    i = (long)(i / sizeof(CALCTYPE));
    M = (long)(i / PROFILEWIDTH);
    /* read in profile */ 
    R = (PROFILECOLUMN *)Malloc( (M+1) * sizeof(PROFILECOLUMN) ); 
    fread( &R[1], sizeof( CALCTYPE ), i, fh);

    fclose(fh);
    *A = R;

    return M;
}

/*--------------------------------------------------------------------*/
/* ma starts at (0;0) */
void ReadCountsFromMA( filename, counts )
    char *filename;
    COUNTCOLUMN *counts;
{
    FILE *fh;
    int i = 0;
    int j, c;
    char *buf;
    SEQTYPE *decode = GetDecode();
    char * buffer = (char *)Malloc( (MAX_SEQUENCE_LENGTH+2) * sizeof(char));

    OpenFile( &fh, filename, "r", "align");
    i = 0;
    while (!feof(fh)) {
      buffer[1] = '\n';
      fgets( &buffer[1], MAX_SEQUENCE_LENGTH, fh); /* read in at position 1; */
      for (j = 1; buffer[j] >= 'A' && buffer[j] <= 'z'; ++j) {
	c = decode[buffer[j]];
	if ( c < PROFILEWIDTH) 
	  counts[j][c] ++;
      }
      ++i;
    }
    fclose( fh );
    free(buffer);
}

/* reads in a MA and puts results into counts */
void ReadMA( filename, width, length, A )
    char *filename;
    int width, length;
    SEQTYPE **A;
{
    FILE *fh;
    int i = 0;
    char buf[NMAX+1];
    char *R;
    SEQTYPE *AA;
    SEQTYPE *decode = GetDecode();
    
    OpenFile( &fh, filename, "r", "align");
    
    R  = (char *)Malloc( width * (length + 1) * sizeof (char));        /* Allocate memory */
    i = 0;
    while (!feof(fh) && (i < width)) {
	fgets( buf, NMAX, fh);
	memcpy( &R[i * length], buf, length);
	i++;
    }
    fclose( fh );

    AA = (SEQTYPE*)Malloc( width * (length + 1) * sizeof (SEQTYPE));        /* Allocate memory */
    for (i = 0; i < (width * length); i++)                                  /* encode */
	AA[i] = decode[R[i]];

    *A = AA;
}

/*--------------------------------------------------------------------*/
int ReadCounts( filename, D )
    char *filename;
    COUNTCOLUMN **D;
{
    FILE *fh;
    int M;
    long i;
    
    COUNTCOLUMN *R;
    
    OpenFile( &fh, filename, "r", "align");
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of count matrix */
    i = (long)(i / sizeof(COUNTTYPE));
    M = (long)(i / PROFILEWIDTH);
    /* read in profile */ 
    R = (COUNTCOLUMN *)Malloc( ((M+1) * PROFILEWIDTH * sizeof (COUNTTYPE)) ); 
    fread( &R[1], sizeof( COUNTTYPE ), i, fh);
    
    fclose(fh);
    *D = R;

    return M;
}

/*--------------------------------------------------------------------*/
/** read a sequence as a string 
    and convert to encoded format.
    The sequence starts at position 1.
 */
int ReadSequence( filename, BB  )
    char *filename;
    SEQTYPE **BB;
{
    FILE *fh;
    int N, i;
    SEQTYPE * decode = GetDecode();
    SEQTYPE * B;

    char * buffer = (char *)Malloc( (MAX_SEQUENCE_LENGTH+2) * sizeof(char));

    buffer[0] = 'X';
    OpenFile( &fh, filename, "r", "align");
    fgets(&buffer[1], MAX_SEQUENCE_LENGTH, fh);
    fclose(fh);

    /* length of sequence (do not count buffer[0] and '\0') */
    N = strlen(buffer) - 1;
    
    /* sometimes there is a newline at the end, discount that */
    while (buffer[N] < 'A' || buffer[N] > 'z' ) N--;

    /* allocate memory and convert sequence */ 
    B = (SEQTYPE*)Malloc( (N + 1) * sizeof (SEQTYPE));           
    for (i = 1; i <= N; i++) 
      B[i] = decode[buffer[i]];

    *BB = B;
    free(buffer);
    return N;
}

/*--------------------------------------------------------------------*/
char * Sequence2String( sequence, lsequence) 
SEQTYPE *sequence;
int lsequence;
{
    char * newstring = Malloc( sizeof(char) * (lsequence + 1));
    int i;

    for (i = 0; i < lsequence; i++) 
	newstring[i] = Encode(sequence[i+1]);

    newstring[lsequence] = '\0';
    return (newstring);
}
/*--------------------------------------------------------------------*/
void WriteFastaFile( filename, title, sequence, lsequence)
char * filename;
char * title;
SEQTYPE * sequence;
int lsequence;
{
    FILE *fh;
    char * newseq = Sequence2String( sequence, lsequence);
    
    OpenFile( &fh, filename, "w", "WriteFastaFile");
    
    fprintf( fh, ">%s\n%s\n", title, newseq);
    
    free( newseq );
    FClose( fh );
}
/*--------------------------------------------------------------------*/
int ReadRealSequence( filename, B  )
    char *filename;
    char **B;
{
    FILE *fh;
    int N;
    long i;

    char *R;

    OpenFile( &fh, filename, "r", "align");
    
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of sequence */
    i = (long)(i / sizeof(char)) - 1; /* subtract eof */
    N = i;
    /* read sequence */ 
    R = (char *)Malloc( ((i+1) * sizeof (char)));
    fread( &R[1], sizeof( char ), i, fh);
    fclose(fh);
    
    *B = R;

    return N;
}
/*--------------------------------------------------------------------*/
int ReadSequenceBinary( filename, B  )
    char *filename;
    SEQTYPE **B;
{
    FILE *fh;
    int N;
    long i;


    SEQTYPE *R;

    OpenFile( &fh, filename, "r", "align");
    
    /* get file size */
    fseek (fh, 0, SEEK_END);
    i = ftell(fh);
    rewind (fh);  
    /* calculate length of sequence */
    i = (long)(i / sizeof(SEQTYPE)) - 1; 
    N = i;
    /* read in profile */ 
    R = (SEQTYPE *)Malloc( ((i+1) * sizeof (SEQTYPE)));
    fread( &R[1], sizeof( SEQTYPE ), i, fh);
    fclose(fh);

    *B = R;
    return N;
}
/*--------------------------------------------------------------------*/

int ReadMatrix( filename, C )
    char *filename;
    MATRIXCOLUMN **C;
{
    FILE *fh;

    MATRIXCOLUMN *R;

    OpenFile( &fh, filename, "r", "align");
    
    /* read in profile */ 
    R = (MATRIXCOLUMN *)Malloc( MATRIXWIDTH * sizeof (MATRIXCOLUMN));
    fread( &R, sizeof( MATRIXCOLUMN ), MATRIXWIDTH, fh);
    fclose(fh);

    *C = R;
    return MATRIXWIDTH;
}

/*--------------------------------------------------------------------*/
void PrintMatrix( matrix )
    MATRIXCOLUMN *matrix;
{
    int i, j;

    for (i = 0; i < MATRIXWIDTH; i++) {
	for (j = 0; j < MATRIXWIDTH; j++) {
	    printf("%6.2f ", matrix[i][j]);
	}
	printf ("\n");
    }
}


/*--------------------------------------------------------------------*/
void PrintAlignment ( ali ) 
ALI ali;
{
    int i;

    for (i = 0; i <= ali.lastindex; i++) {
	printf("%5i%5i%10.2f\n", 
	       ali.align_i[i], ali.align_j[i], ali.align_s[i]);
    }
}

/*--------------------------------------------------------------------*/
void PrintProfile ( A, M )
PROFILECOLUMN *A;
int M;
{
    int i,j;
    
    for (i = 1; i <= M; i++)
      {
	printf("%2i\t", i);
	for (j = 0; j < PROFILEWIDTH; j++)
	  {
	    printf("%6.2f", A[i][j]);
	  }
	printf("\n");
    }
}

/*--------------------------------------------------------------------*/
void PrintFrequencies ( A, M )
FREQUENCYCOLUMN *A;
int M;
{
    int i,j;
    
    for (i = 1; i <= M; i++) {
	printf("%2i\t", i);
	for (j = 0; j < PROFILEWIDTH; j++) {
	    printf("%6.2f", A[i][j]);
	}
	printf("\n");
    }
}

/*--------------------------------------------------------------------*/
void PrintCounts ( A, M )
COUNTCOLUMN *A;
int M;
{
    int i,j;
    
    for (i = 1; i <= M; i++) {
	printf("%2i\t", i);
	for (j = 0; j < PROFILEWIDTH; j++) {
	    printf("%6i", A[i][j]);
	}
	printf("\n");
    }
}

/*--------------------------------------------------------------------*/
void PrintSequence (B, N) 
SEQTYPE *B;
int N;
{
    int i;
    
    for (i = 1; i <= N; i++) {
	printf("%4i", B[i]);
    }
    printf("\n");
}    

/*--------------------------------------------------------------------*/
void PrintMA (B, M, N) 
SEQTYPE *B;
int M,N;
{
    int i,j;
    
    for (i = 0; i < M; i++) {
	for (j = 0; j < N; j++) {
	    printf("%4i", B[j + i * N]);
	}
	printf("\n");
    }
}    

/*--------------------------------------------------------------------*/
void PrintDots( row, col, score, M )
int *row, *col;
CALCTYPE *score;
int M;
{
    int i;
    for (i = 0; i < M; i++) {
	printf("%i %i %i %5.2f\n", i, row[i], col[i], score[i]);
    }
}

/*--------------------------------------------------------------------*/
void WriteDots( row, col, score, M, filename )
int *row, *col;
CALCTYPE *score;
int M;
char * filename;
{
    FILE * fh;
    OpenFile( &fh, filename, "w", "WriteProfile");
    fwrite( row, sizeof(int), M, fh);
    fwrite( col, sizeof(int), M, fh);
    fwrite( score, sizeof(CALCTYPE), M, fh);
    FClose( fh );
}

/*--------------------------------------------------------------------*/
void WriteProfile ( A, M, filename )
PROFILECOLUMN *A;
int M;
char *filename;
{
    FILE * fh;
    
    OpenFile( &fh, filename, "w", "WriteProfile");
    fwrite( A[1], sizeof(PROFILECOLUMN), M, fh);
    FClose( fh );
}

/*--------------------------------------------------------------------*/
void WriteCounts ( A, M, filename )
COUNTCOLUMN *A;
int M;
char *filename;
{
    FILE *fh;
    
    OpenFile( &fh, filename, "w", "WriteCounts");
    fwrite( A[1], sizeof(COUNTCOLUMN), M, fh);
    FClose( fh );

}

