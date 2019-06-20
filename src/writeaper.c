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

void writeaper(double *pix, int W, int H, char *img)
{

  fitsfile *thptr;
  int status = 0; /*  CFITSIO status value MUST be initialized to zero! */
  int i;
  long firstpix[2] = {1,1};
  double *Row, *ptr;
  int bitpix=-64;
  int naxis=2;
  long *naxes;

  naxes = (long *) malloc(2 * sizeof(long));
  naxes[0]=W;
  naxes[1]=H;

  /* Allocate memory for storing double image row */
  if (fits_create_file(&thptr, "!aper.fits", &status)) {
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
  }
  
  if (fits_create_img(thptr, bitpix, naxis, naxes, &status)) {
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
  }
  if(fits_write_pix(thptr, TDOUBLE, firstpix, W*H, pix, &status)) {
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
  }          
  
  if(fits_close_file(thptr, &status)) {
    fits_report_error(stderr, status);
    exit(EXIT_FAILURE);
  }  

  free(naxes);

}

