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
//#define DEBUG

float a,b,theta;

void moments(int Lx, int Ly, int xmin, int ymin, float x0, float y0, double *pixels, double *pixels_rms, int *pixels_seg, int *pixels_flag, float bg, int id, float pa)
{

  double xc,yc,xcc,ycc,d2;
  int i,j,n,dl;
  double x2m,y2m,xym,s,r2,f;

  int clipmax=100;
  int clip, TOL=5;
  float fac_sigma_mom=2.0; // Level of sigma-clipping

  xc=x0-xmin;
  yc=y0-ymin;
  int ic=(int)xc;
  int jc=(int)yc;
  r2=MAX(xc,yc);
  r2=r2*r2;

  // compute second moments
  
  int usesegm=0;
  for (j=jc-1;j<=jc+1;j++)
    for (i=ic-1;i<=ic+1;i++) {
      if (pixels_seg[j*Lx+i]==id)
	{
	  usesegm=1;
#ifdef DEBUG
      printf("Using segmentation\n");
#endif
	  break;
	}
    }
  
  if (usesegm==1) 
    {
      // Just use segmented pixels
      xcc=0.0; ycc=0.0;
      x2m=0.0;
      y2m=0.0;
      xym=0.0;
      s=0.0;
      n=0;

      for (j=0;j<Ly;j++){
	for (i=0;i<Lx;i++){

	  if (pixels_seg[j*Lx+i]==id)
	    if (pixels_flag[j*Lx+i]==0){
	      {
		if ((i-xc)*(i-xc)+(j-yc)*(j-yc)<=r2)
		  {
		    f=pixels[j*Lx+i]-bg;
		    xcc+=f*i;
		    ycc+=f*j;
		    x2m+=f*i*i;
		    y2m+=f*j*j;
		    xym+=f*i*j; 
		    s+=f;
		    n+=1;
		  }
	      }
	    }
	}
      }
    }
  else // WARNING - THIS SECTION IS TENTATIVE, CURRENTLY NO GOOD RESULTS
    {
      // No segmentation; select and use only pixels above fac_sigma_mom
      clip=0;
      s=0.0;
      dl=4;
      while (clip>=0)
	{
	  xcc=0.0; ycc=0.0;
	  x2m=0.0;
	  y2m=0.0;
	  xym=0.0;
	  s=0.0;	  
	  n=0;

	  //for (j=2;j<=Ly-2;j++){
	  //  for (i=2;i<=Lx-2;i++){
	  for (j=yc-dl;j<=yc+dl;j++){
	    for (i=xc-dl;i<=xc+dl;i++){
	      if (pixels_flag[j*Lx+i]==0){
		d2=(i-xc)*(i-xc)+(j-yc)*(j-yc);
		if (d2<=r2)
		  {
		    f=pixels[j*Lx+i]-bg;
		    int ii, jj;
		    double fff=0.0, fffrms=0.0;
		    for (ii=i-2;ii<=i+2;ii++)
		      for (jj=j-2;jj<=j+2;jj++)
			{
			  fff+=pixels[jj*Lx+ii]-bg;
			  fffrms+=pixels_rms[jj*Lx+ii];
			}
		    if (fff>fac_sigma_mom*fffrms) 
		      {
			xcc+=f*i;
			ycc+=f*j;
			x2m+=f*i*i;
			y2m+=f*j*j;
			xym+=f*i*j;
			s+=f;
			n+=1;
		      }
		  }
	      }
	    }
	  }

	  xcc=(x0-xmin)*s; ycc=(y0-ymin)*s;
	  if (n>TOL)
	    {clip=-1;}
	  else if (clip==clipmax)
	    {
#ifdef VERBOSE
	      printf("Reached clipmax %d / %d\n",clip,clipmax);
#endif	      
	      clip=-1;
	    }
	  else
	    {
#ifdef VERBOSE
	      printf("Moments - too few pixels (%d); fac_sigma_mom -> %f\n",n,fac_sigma_mom-0.1);
#endif
	      fac_sigma_mom-=0.1;
	      clip+=1;
	    }	    
	}
    }
  
  // baricenter
  xcc=xcc/s;
  ycc=ycc/s;

#ifdef DEBUG
  printf("^^^ %d %lf %lf %lf %lf %lf -> %lf\n",n,x2m,s,x2m/s,xcc,xcc*xcc,x2m/s-xcc*xcc);
#endif  
  // second moments
  x2m=x2m/s-xcc*xcc;
  y2m=y2m/s-ycc*ycc;
  xym=xym/s-xcc*ycc;

#ifdef DEBUG
  printf("baricenter %lf %lf (%lf %lf)\n",xcc+xmin,ycc+ymin,x0,y0);
  printf("moments %lf %lf %lf\n",x2m,y2m,xym);
#endif  

  a=MAX(0.25,0.5*(x2m+y2m)+sqrt(pow(0.5*(x2m-y2m),2)+xym*xym));
  b=MAX(0.01*a,0.5*(x2m+y2m)-sqrt(pow(0.5*(x2m-y2m),2)+xym*xym));
  a=sqrt(a);
  b=sqrt(b);

  if (pa>=500.0)
    {
      theta=0.5*atan(2.0*xym/(x2m-y2m+1.e-8));
      if (xym*theta<0.0) {theta-=HALFPI;}
      if (theta<-HALFPI) {theta+=2.0*HALFPI;}
      if (theta>HALFPI) {theta-=2.0*HALFPI;}
    }
  else
    {
      theta=pa/57.2958;
    }
    
#ifdef VERBOSE
  printf("ooo> a %f b %f theta %f\n",a,b,57.2958*theta);
#endif
  
}
