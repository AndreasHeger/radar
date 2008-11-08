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
$Id: align_lib.h,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository


*/
#define __ALIGN_LIB_H

/* ------------ alignment functions ----------- */
void lps_repeat();
void lssb_align();
void lsp_repeat();
void gss_align();
void gpp_align();
void gps_align();
void lpp_align();
void lps_align();
CALCTYPE lps_align_score();
void lss_align();
void lsp_repeat_combi();
void lsp_repeat_insert();
void dotalign();
void dotalign_tube();
void dotalign_gtube();
void dotalign_gonnet();
int  getmaxdiag();
int  getnbestdiag();
int  getnbestdiagsum();
int  getnbestdiagpostrace();
void mask_dots();
int  Fasta2Dots();
int  Fasta2DotsCorrect();
int  Fasta2Dotfile();
int  FastaFiles2DotsCorrect();
int  FastaFiles2DotsCorrectSaveSpace();

void score_dots_profile();
void add_diagonal();

void Counts2Profile();
void Counts2Frequencies();
void MA2Profile();
void MA2Counts();
double MA_Analyse();

