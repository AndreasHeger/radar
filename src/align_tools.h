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
$Id: align_tools.h,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#define __ALIGN_TOOLS_H

/* ----------- macros -----------------------*/

#define gap(k)   ((k) <= 0 ? 0 : g +h *(k))	/* k-symbol indel cost */
#define gapi(k)  ((k) <= 0 ? 0 : gi+hi*(k))	/* k-symbol indel cost */
#define gapj(k)  ((k) <= 0 ? 0 : gj+hj*(k))	/* k-symbol indel cost */


#define DEL(k)				\
{ if (last < 0)				\
    last = sapp[-1] -= (k);		\
  else					\
    last = *sapp++ = -(k);		\
}
						/* Append "Insert k" op */
#define INS(k)				\
{ if (last < 0)				\
    { sapp[-1] = (k); *sapp++ = last; }	\
  else					\
    last = *sapp++ = (k);		\
}

#define REP { last = *sapp++ = 0; }		/* Append "Replace" op */

typedef struct {
    int i;
    int j;
} CELL;

#define NCOMPONENTS  9
#define SCALEFACTOR  10.0

void ProfileInitialize();
void RegularizeColumn();
void RegularizeColumnProbability();
SEQTYPE *GetDecode();
SEQTYPE Decode();
char Encode();

void Counts2Profile();
void Counts2Frequencies();
void MA2Profile();
void MA2Counts();

/* ------------ helper functions --------------- */
CALCTYPE DotProduct();
CALCTYPE ProfileScore();
void TraceBack();
void TraceBackDotalignment();
void ConvertAli();
int  WasGapHere();

#define INDEX(i,j) (i*O+j)
#define STACKEMPTY 0    
#define NOMATCH -999
#define INSERTION -1
#define DELETION  -2
#define MATCH     -3
#define ROLL      -4

