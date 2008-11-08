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
$Id: toolbox.c,v 1.1.1.1 2003/10/01 07:29:24 aheger Exp $

ChangeLog:

15.3.2001	heger	Added to repository

*/
#include <stdlib.h>
#include <stdio.h>
#include "toolbox.h"
#include "align.h"

void memerror()
{
printf("Error allocating memory\n");
exit(-1);
}

MALLOCRETURN *mymalloc(x)
long x;
{
    MALLOCRETURN *mem;

#ifdef DEBUG 
    printf ("Allocating %i bytes of memory\n", x);
#endif
    mem = (MALLOCRETURN *)malloc(x);
	
    if (!mem)
	memerror();
    
    return (MALLOCRETURN *)mem;
} 

void uppercase(ch) 
    char *ch; { /* make ch upper-case */    
    *ch = (islower(*ch) ?  toupper(*ch) : (*ch)); 
}  /* uppercase */                            

/*--------------------------------------------------------------------*/

void OpenFile(fp,filename,mode,application)
    FILE **fp;
char *filename;
char *mode;
char *application;
{
  FILE *of;
  char file[100];
  strcpy(file,filename);
  while (1){
    of = fopen(file,mode);
    if (of)
      break;
    else {
	printf("Error: file %s not found\n", file);
	exit(0);
	switch (*mode){
      case 'r':
        printf("%s:  can't read file %s\n",application,file);
        file[0]='\0';
        while (file[0] =='\0') {
            printf("Please enter a new filename>");
	    gets(file); }
        break;
      case 'w':
        printf("%s: can't write %s\n",application,file);
        file[0] = '\0';
        while (file[0] =='\0'){
            printf("Please enter a new filename>");
            gets(file);}
        break;
      }
    }
  }
  *fp=of;
}

/* x is sorted by y */
 
void QuickSortIndicesDouble (from, to, x, y)
int from, to;
int      *x; 
CALCTYPE *y;

{
    int lastsmall, t, c1, i;
    int m; /* value of median */
    CALCTYPE c;
    

/*     printf(" from %i, to %i \n", from, to); */
/*     for (i = 0; i < 12; i++) printf("%5.2f - ", y[x[i]]); */
/*     printf("\n"); */

    if (from < to ) {
	c1 = (from + to) / 2;         
	t = x[from]; x[from] = x[c1]; x[c1] = t;
	c = y[x[from]]; /* choose median value to compare to */
	lastsmall = from;
	for (i = from + 1; i <= to; i++) {
	    if ( y[x[i]] < c) { /* swap lastsmall and i */
		lastsmall++;
		t = x[lastsmall]; x[lastsmall] = x[i]; x[i] = t;
	    }
	}
	t = x[from]; x[from] = x[lastsmall]; x[lastsmall] = t; 
	
	m = lastsmall;
	QuickSortIndicesDouble( from, m, x, y);
	QuickSortIndicesDouble( m+1, to, x, y);
    }
}
