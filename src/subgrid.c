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
//#define INTERPSUBGRID

float subgrid(double *fff, double **dw, int ii, int jj, float xm, float ym, float x0, float y0, int binsubpix, double dsubpix, double* pixels, int Lx, int Ly, double R2, float aa, float bb, double cosPPAA, double sinPPAA)
{

  double sfacarea=0.0;
  double fdw;
  double ffftot=0.0;
  int ni,nj;
  double dx, dy, dsx, dsy;
  int iii, jjj, i, j;

  dx=xm-x0;
  dy=ym-y0;

  // First, interpolate to obtain subgrid fluxes
  for (iii=0; iii<binsubpix*binsubpix; iii++) fff[iii]=0.0;

  for (jjj=0;jjj<binsubpix;jjj++)
    for (iii=0;iii<binsubpix;iii++)
      {
#ifdef INTERPSUBGRID
	fdw=0.0;
	for (nj=0;nj<3;nj++)
	  for (ni=0;ni<3;ni++)
	    {
	      if ((nj+jj-1>=0)&&(ni+ii-1>=0)&&(nj+jj-1<Ly)&&(ni+ii-1<Lx))
		{
#ifdef VERBOSE
		  printf("+> %d %d %d %d %d %d\n",ni,nj,ni+ii-1,nj+jj-1,Lx,Ly);
#endif	      
		  // compute weight as the inverse of distance from the center of this pixel
		  fdw+=dw[nj*3+ni][jjj*binsubpix+iii];
		  //printf("___ %d %d %d %d %f\n",iii,jjj,ni,nj,dw[nj*3+ni][jjj*binsubpix+iii]);
		  fff[jjj*binsubpix+iii]+=dw[nj*3+ni][jjj*binsubpix+iii]*pixels[(nj+jj-1)*Lx+(ni+ii-1)];
		}
	    }
	fff[jjj*binsubpix+iii]/=fdw;
#else
        fff[jjj*binsubpix+iii]=pixels[jj*Lx+ii]/(float)(binsubpix*binsubpix);
#endif
	ffftot+=fff[jjj*binsubpix+iii];
      }

#ifdef INTERPSUBGRID
#ifdef DEBUG
  printf("[------------]\n");
  for (jjj=0;jjj<binsubpix;jjj++)
    {
      for (iii=0;iii<binsubpix;iii++)
        {
          printf("%f ",fff[jjj*binsubpix+iii]);
        }
      printf("\n");
    }
  printf("[------------]\n");
#endif
#endif

  // Then, only include subpixels within radius
  if (aa<0.0)
    {
      for (jjj=0;jjj<binsubpix;jjj++)
	for (iii=0;iii<binsubpix;iii++)
	  {
	    dsx=dx+((float)iii+0.5)*dsubpix;
	    dsy=dy+((float)jjj+0.5)*dsubpix;
	    if (dsx*dsx+dsy*dsy<=R2) {
	      sfacarea+=fff[jjj*binsubpix+iii]/ffftot;
	    }
	  }
    }
  else
    {
      for (jjj=0;jjj<binsubpix;jjj++)
	for (iii=0;iii<binsubpix;iii++)
	  {
	    
	    dsx=(dx+((float)iii+0.5)*dsubpix)*cosPPAA+(dy+((float)jjj+0.5)*dsubpix)*sinPPAA;
	    dsy=-(dx+((float)iii+0.5)*dsubpix)*sinPPAA+(dy+((float)jjj+0.5)*dsubpix)*cosPPAA;
	    
	    if ((dsx*dsx)/(aa*aa)+(dsy*dsy)/(bb*bb)<=1.0) {
	      sfacarea+=fff[jjj*binsubpix+iii]/ffftot;
#ifdef VERBOSE
	      //printf("* %f %f %f\n",sfacarea,fff[jjj*binsubpix+iii],ffftot);    
#endif
	    }
	  }
    }


#ifdef VERBOSE
  printf("... %f %f\n",sfacarea,pixels[nj*Lx+ni]);
#endif

  return(sfacarea);
}


//---------------------------------------------------------------------------//

