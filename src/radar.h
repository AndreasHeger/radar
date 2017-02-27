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
$Id: radar.h,v 1.2 2005/02/17 23:13:32 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/

#ifndef _RADAR_H
#define _RADAR_H

#if HAVE_CONFIG_H
#include <config.h>
#endif


#include "align.h"
#include "align_tools.h"
#include "align_io.h"
#include "toolbox.h"
#include "align_lib.h"

#define MAXNREPEATS 500          /* you won't believe it, but 100 might be too small */
#define UPPER 1
#define INFINITE 1000
#define RINFINITE 1000.0
#define MAX_NITERATIONS 100
#define MAXRESULTSIZE 1000000


typedef struct {
    int * row;
    int * col;
    CALCTYPE *score ;
    int ndots;
    int ndotsdiagonal;
} DOTS;

typedef struct {
    int from, to;
    CALCTYPE score;
} SEQUENCEWINDOW;

typedef struct {
    int* ali;
    CALCTYPE *scores;
    CALCTYPE score;
    int length;
} ARRAYALI;
    

typedef enum {
    LL0, LL1, LL2, LL3, LL4, LL5
} VERBOSITY;

void DumpAlignment();
void DumpArrayAlignment();
void AddDiagonal ();
ALI *AllocateAliMemory ();
ARRAYALI *AllocateArrayAliMemory ();
void FreeArrayAliMemory();
COUNTCOLUMN *AllocateCountsMemory ();
PROFILECOLUMN *AllocateProfileMemory ();
void FreeCountsMemory () ;
void FreeProfileMemory ();
void FreeAliMemory () ;
void FreeDotsMemory () ;
DOTS* AllocateDotsMemory ();
DOTS *CopyDots();
ALI *CopyAli();
ARRAYALI *CopyArrayAli();
DOTS* MaskDots ();
int FoldAlignment ();
ARRAYALI * Ali2ArrayAli ();
SEQTYPE  *CopySequence();
SEQUENCEWINDOW FindBestWindow ();
SEQUENCEWINDOW ReduceDiagonal ();
SEQUENCEWINDOW GetLastRepeat ();
SEQUENCEWINDOW GetFirstRepeat ();
void PrintPrettyMA ();
CALCTYPE ShiftAli ();
CALCTYPE * GetGapVote ();
int VoteGap ();
int VotePushGap ();
CALCTYPE PushGapVote (); 
CALCTYPE UltimateRegister (); 
CALCTYPE NewRegister (); 
void NewestRegister (); 
int CountAndScoreRepeats();
void CalculateZScore();


int radar_run_from_files(const char *, const char *, const char *, const char *, unsigned int);
void radar_setLogLevel(int);

#endif
