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
$Id: counts2profile.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "align.h"


/*-------------------------------*/

void 
Counts2Profile (A, B, M )
COUNTCOLUMN   *A;
PROFILECOLUMN *B;
int            M;
{
    int       i;

#ifdef DEBUG
    printf("----------- Counts2Profile called ---------- \n");
    PrintCounts( A, M);
#endif
    /* Initialize some variables */
    ProfileInitialize();
    
    /* profile starts at 1, not at 0 !!!*/
    for (i = 1; i <= M; i ++) {
	RegularizeColumn( A[i], B[i]);
    }
}	

