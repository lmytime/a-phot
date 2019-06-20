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

//#define VERBOSE

#include <string.h>
#include <stdio.h>
#include "fitsio.h"
#include <aphot.h>
#include "globals.h"

int phot_flag(float bad_pix, float aper_area, float contamin_pix, double ffff, int blend_obj, int pix_max, float x0, float y0, double RRR, double Rmax, double cosPA, double sinPA, int W, int H)
{

#ifdef VERBOSE
  printf(">>> %f %f %f %lf %d %d %f %f %lf %lf %lf %lf %d %d\n",bad_pix,aper_area,contamin_pix,ffff,blend_obj,pix_max,x0,y0,RRR,Rmax,cosPA,sinPA,W,H);
#endif						
  
  int obj_flag=0;
  
  // Check object flagging
  if ((bad_pix/aper_area >= 0.1) || (contamin_pix>ffff))
    {
      // contamination / bad pixels
      obj_flag=+1;
    }
  
  if (blend_obj>0)
    {
      // blending
      obj_flag+=2;
    }
  
  if (pix_max>2)
    {
      // saturation
      obj_flag+=4;
    }
  
  // border
  if ((x0+RRR*cosPA>=(float)W) || (x0-RRR*cosPA<=0.0) || (y0+RRR*sinPA>=(float)H) || (y0-RRR*sinPA<=0.0))
    {obj_flag+=8;}
  else
    {
      if ((x0+Rmax>=(float)W) || (x0-Rmax<=0.0) || (y0+Rmax>=(float)H) || (y0-Rmax<=0.0))
	{obj_flag+=8;}
    }

  return(obj_flag);
}
