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

int xmin,xmax,ymin,ymax,Lx,Ly,Area, W,H;
long fpixel[2], lpixel[2], inc[2] = {1,1};
fitsfile *fptr, *fptrrms, *fptrseg, *fptrflag;
double *pixels, *pixels_orig, *pixels_rms;
int *pixels_seg, *pixels_seg_orig, *pixels_flag;

void readpixels(int xc, int yc, int rrr, int useseg, int useflag)
{
  int status=0;
  int i;
  
  // Define region of interest
  xmin=MAX(0,xc-rrr-1);
  xmax=MIN(xc+rrr+1,W-1);
  ymin=MAX(0,yc-rrr-1);
  ymax=MIN(yc+rrr+1,H-1);
  Lx=(xmax-xmin+1);
  Ly=(ymax-ymin+1);
  Area=Lx*Ly;
  
  // Read in pixel values
  fpixel[0]=xmin+1; fpixel[1]=ymin+1;
  lpixel[0]=xmax+1; lpixel[1]=ymax+1;

  
  // scientific image
  pixels = (double *) malloc(Area * sizeof(double));
  fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, NULL,
		   pixels, NULL, &status);

  // rms image
  pixels_rms = (double *) malloc(Area * sizeof(double));
  fits_read_subset(fptrrms, TDOUBLE, fpixel, lpixel, inc, NULL,
		   pixels_rms, NULL, &status);

  // segmentation
  pixels_seg = (int *) malloc(Area * sizeof(int));
  if (useseg==1)
    {
      fits_read_subset(fptrseg, TINT, fpixel, lpixel, inc, NULL,
		       pixels_seg, NULL, &status);
    }
  else
    {
      for (i=0;i<Area;i++) pixels_seg[i]=0;
    }

  //flag image
  pixels_flag = (int *) malloc(Area * sizeof(int));
  if (useflag==1)
    {
      fits_read_subset(fptrflag, TINT, fpixel, lpixel, inc, NULL,
		       pixels_flag, NULL, &status);
    }
  else
    {
      for (i=0;i<Area;i++) pixels_flag[i]=0;
    }


  // Now again, for pixel replacing
  
  // scientific image
  pixels_orig = (double *) malloc(Area * sizeof(double));
  fits_read_subset(fptr, TDOUBLE, fpixel, lpixel, inc, NULL,
		   pixels_orig, NULL, &status);
  
  // segmentation
  pixels_seg_orig = (int *) malloc(Area * sizeof(int));
  if (useseg==1)
    {
      fits_read_subset(fptrseg, TINT, fpixel, lpixel, inc, NULL,
		       pixels_seg_orig, NULL, &status);
    }
  else
    {
      for (i=0;i<Area;i++) pixels_seg_orig[i]=0;
    }
  
}


//---------------------------------------//

void freepixels()
{
  free(pixels_flag);  
  free(pixels_seg);
  free(pixels_rms);
  free(pixels);
 
  free(pixels_seg_orig);
  free(pixels_orig);

}
