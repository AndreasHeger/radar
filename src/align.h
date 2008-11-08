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
$Id: align.h,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository

*/

#ifndef __ALIGN_H
#define __ALIGN_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ------------ memory -------------------- */
#define NMAX 30000              /* maximum sequence length   */
#define MAXALI_SIZE 30000       /* maximum size of alignment */
#define MAXNDOTS 2000000        /* maximum number of dots (2000000 for titin, use less in other cases) */
#define MAXLENBUFFER 1000       /* maximum size of one line */
/* ------------ miscallenous -------------- */

#define MASKCODE      20	           /* character corresponding to a masked residues */
#define UNKNOWN_CODE 'X'             /* character corresponding to a residue with unknown code */
#define GAPCODE       21             /* code corresponding to gap */

#define PROFILEWIDTH 20              /* width of profile, 4 for NA, 20 for AA */
#define CALCTYPE     double          /* type of profile and dots*/
#define SEQTYPE      unsigned char   /* type of sequence integers */
#define COUNTTYPE    unsigned int    /* type of counts */
#define MATYPE       unsigned char   /* type of multiple alignment */
#define MATRIXWIDTH  24              /* width of substitution matrices */       

/* --------------- typedefs ------------------*/


typedef CALCTYPE  PROFILECOLUMN[PROFILEWIDTH];
typedef double    FREQUENCYCOLUMN[PROFILEWIDTH];
typedef COUNTTYPE COUNTCOLUMN[PROFILEWIDTH];
typedef CALCTYPE  MATRIXCOLUMN[MATRIXWIDTH];

typedef struct {
    int length;                        /* length of alignment      */
    int lastindex;                     /* last index in ali-arrays */
    int *align_i, *align_j;            /* aligned residues         */
    CALCTYPE *align_s;                 /* score array              */
    CALCTYPE score;                    /* alignment score          */
} ALI;   

typedef int* INTP;

#endif




