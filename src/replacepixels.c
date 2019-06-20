/*
  Copyright 2016 Emiliano Merlin

  This file is part of A-PHOT.

  A-PHOT is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  A-PHOT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with A-PHOT.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <aphot.h>
#include "globals.h"

//#define VERBOSE

double *pixels, *pixels_rms;
int *pixels_seg, *pixels_flag;

void replace_sym(int xxx, int yyy, int id, int Lx, int Ly, double RMSTOL)
{
  int i,j,idx,xn,yn;
  float fff;

#ifdef VERBOSE
  printf(" pixreplacing >>>%d\n",id);
#endif

  for(j=0;j<Ly;j++)
    {
    for(i=0;i<Lx;i++)
      {
	if ((pixels_seg[j*Lx+i]!=0 && pixels_seg[j*Lx+i]!=id) || pixels_flag[j*Lx+i]!=0 || pixels_rms[j*Lx+i]>RMSTOL)
	  {

	    xn=(2*xxx-i);
	    yn=(2*yyy-j);

	    //Check if it's outside the box
	    if ((xn>=Lx && yn>=Ly) || (xn<0 && yn<0))
	      {
		idx=j*Lx+i; // no symmetry possible, it will be flagged
		fff=0.0;
	      }
	    else if ((xn>=Lx && yn<Ly) || (xn<0 && yn>=0))
	      {
		idx=yn*Lx+i;
		fff=pixels[idx];
	      }
	    else if ((yn>=Ly && xn<Lx) || (yn<0 && xn>=0))
	      {
		idx=j*Lx+xn;
		fff=pixels[idx];
	      }
	    else
	      {
		idx=yn*Lx+xn; // symmetric index
		fff=pixels[idx];
	      }
	
	    if ((pixels_seg[idx]==0 || pixels_seg[idx]==id) && pixels_flag[idx]==0 && pixels_rms[idx]<=RMSTOL)
	      {
#ifdef VERBOSE
		printf("      >> replacing pix %d: %d %d %d %f -> %f\n",j*Lx+i,idx,pixels_seg[j*Lx+i],pixels_flag[j*Lx+i],pixels[j*Lx+i],fff);
#endif
		pixels[j*Lx+i]=fff; // symmetric substitution
		pixels_rms[j*Lx+i]=pixels_rms[idx]; // symmetric substitution
		pixels_flag[j*Lx+i]=0;
		pixels_seg[j*Lx+i]=0;
	      }
#ifdef VERBOSE
	    else
	      {
		printf("NOT REPLACING! %d %d\n",pixels_seg[idx],pixels_flag[idx]);
	      }
#endif
	    
	  }
      }
    }
}


