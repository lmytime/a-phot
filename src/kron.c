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

//#define RAW
//#define WRITEKRON
#define f_kron 1.0 //0.25
#define r 6.0

double kron(int Lx, int Ly, int xmin, int ymin, float x0, float y0, float a, float b, float amax, double *pixels, int *pixels_seg, int *pixels_flag, float bg, int id, float fmin_kron, float fmax_kron)
{
  double R1=0.0;
  double dR1=0.0;
  float fac2_r=36.0;
  int i,j,flag;
  double x,y,x1,y1,cxx,cyy,cxy,r2,kr2;
  FILE *log1;
  char logname[BUFFERSIZE];
  
  /*
  if (r*sqrt(a*b)>amax)
    {
      cxx=a*a*sinPA*sinPA+b*b*cosPA*cosPA;
      cyy=a*a*cosPA*cosPA+b*b*sinPA*sinPA;
      cxy=2.0*sinPA*cosPA*(b*b-a*a);
      //cxx=0.5*(1.0+cos2PA)/(a*a) + 0.5*(1.0-cos2PA)/(b*b);
      //cyy=0.5*(1.0+cos2PA)/(b*b) + 0.5*(1.0-cos2PA)/(a*a);
      //cxy=2*cosPA*sinPA*(1.0/(a*a)-1.0/(b*b));
      kr2=r*r;
      flag=1;
    }
  else
    {
      cxx=1.0;
      cyy=1.0;
      cxy=0.0;
      kr2=amax*amax;
      flag=2;
    }

#ifdef RAW
  cxx=1.0;
  cyy=1.0;
  cxy=0.0;
  kr2=a*a;
#endif
  
  for (j=0;j<Ly;j++)
    for(i=0;i<Lx;i++) // Loop on pixels of this object
      {

	if ((pixels_seg[j*Lx+i]==id) || (pixels_seg[j*Lx+i]==0)) 
	  {
	    x=(float)(xmin+i);
	    y=(float)(ymin+j);
	    
#ifdef RAW
	    r2=(x-x0)*(x-x0)+(y-y0)*(y-y0);
	    if (r2<=kr2*a*a) //(sqrt(r2)<MAX(Lx,Ly)/2.0) //(r2<=kr2*a*a)
	      {
		R1+=sqrt(r2)*pixels[j*Lx+i];
		dR1+=pixels[j*Lx+i];
		//printf(">>> %lf\n",pixels[j*Lx+i]);
	      }
#else	
	    r2=cxx*(x-x0)*(x-x0)+cyy*(y-y0)*(y-y0)+cxy*(x-x0)*(y-y0);
	    if (r2<=kr2*a*a)
	      {
		//printf("%f %f\n",r2,kr2*a*a);
		R1+=sqrt(r2)*pixels[j*Lx+i];
		dR1+=pixels[j*Lx+i];
	      }
#endif
	    }
      }
	
#ifdef RAW
  R1=2.5*R1/dR1;
#else
  if (dR1>0.0) 
    {
      R1=f_kron*R1/dR1;
    }
  else
    {
      R1=f_kron*a;
    }
#endif
  printf("xxxxooo>>>>> %f\n",R1);
  */
  /*
  // First approach, NOT GOOD
  R1=0.0; dR1=0.0;
  for(i=0;i<Lx;i++){
    cxx=a*a*sinPA*sinPA+b*b*cosPA*cosPA;
    cyy=a*a*cosPA*cosPA+b*b*sinPA*sinPA;
    cxy=2.0*sinPA*cosPA*(b*b-a*a);
    cxx=cxx/(a*a+b*b);cyy=cyy/(a*a+b*b);cxy=cxy/(a*a+b*b);
    j=(int)(Ly/2);
    x=(float)(xmin+i);
    y=(float)(ymin+j);
    r2=cxx*(x-x0)*(x-x0)+cyy*(y-y0)*(y-y0)+cxy*(x-x0)*(y-y0);
    //printf("%f %f\n",x-x0,pixels[j*Lx+i]);
    R1+=r2*pixels[j*Lx+i];
    dR1+=sqrt(r2)*pixels[j*Lx+i];
  }
  R1=2.5*f_kron*R1/dR1;
  printf(">>>>>> %f\n",R1);
  */

  /*
  // New approach, seems to work
  R1=0.0; dR1=0.0;
  printf("..... %d %d\n",Lx,xmin);
  j=(int)(Ly/2);
  for(i=0;i<Lx;i++){
    x=(float)(xmin-1+i);
    y=(float)(ymin-1+j);
    float x1,y1;
    x1=(x-x0)*cosPA+(y-y0)*sinPA;
    y1=-(x-x0)*sinPA+(y-y0)*cosPA;
    r2=x1*x1+y1*y1;
    int ii,jj;
    ii=(int)(Lx/2.0+x1);
    jj=(int)(Ly/2.0-y1);
    if (jj*Lx+ii<Lx*Ly){
      R1+=r2*pixels[jj*Lx+ii];
      dR1+=sqrt(r2)*pixels[jj*Lx+ii];
    }
    //printf("%d %d - %d %d %d %f %f %f %f\n",i,j,ii,jj,jj*Lx+ii,r2,pixels[jj*Lx+ii],R1,dR1);
    printf("%d %d\n",xmin+ii,ymin+jj);
  }
  R1=2.5*f_kron*R1/dR1;
  printf("ooo>>>>>> %f\n",R1);
  */

  
  // Try in two dimensions
  R1=0.0; dR1=0.0;
  //printf("..... %d %d %f - %d\n",Lx,Ly,b/a,id);
  for(j=0;j<Ly;j++){
    for(i=0;i<Lx;i++){
      x=(float)(xmin-1+i);
      y=(float)(ymin-1+j);
      float x1=0.0;
      float y1=0.0;
      x1=x-x0; //(x-x0)*cosPA+(y-y0)*sinPA; //*(a*a/(b*b));
      y1=y-y0; //-(x-x0)*sinPA+(y-y0)*cosPA; //*(a*a/(b*b));
      r2=x1*x1+y1*y1;
      int ii,jj;
      ii=(int)(Lx/2.0+x1);
      jj=(int)(Ly/2.0+y1);
      if (jj*Lx+ii<Lx*Ly){
	if (jj*Lx+ii>=0){
	  if ((pixels_seg[jj*Lx+ii]==id)||(pixels_seg[jj*Lx+ii]==0)){
	    if (pixels_flag[jj*Lx+ii]==0){
	      if (r2<=fac2_r*a*a){
		R1+=sqrt(r2)*(pixels[jj*Lx+ii]-bg);
		dR1+=(pixels[jj*Lx+ii]-bg);
		//printf("%f %f - %d %d %d %f %f %f\n",x1,y1,ii,jj,jj*Lx+ii,pixels[jj*Lx+ii],R1,dR1);
	      }
	    }
	  }
	}
      }
      
    }
  }
  R1=MAX(f_kron*R1/dR1*2.5,fmin_kron*a);
  if (R1>fmax_kron*a)
    {
      R1=fmax_kron*a;
#ifdef VERBOSE
      printf("WARNING: ID %d - cutting Kron radius to 20a=%f\n",id,R1);
#endif      
      //printf("%f %f\n",x0,y0);
    }

#ifdef WRITEKRON
  float s=0.0;
  snprintf(logname,BUFFERSIZE,"aper-%d.txt",id);
  log1 = fopen(logname, "w");
  for (j=0;j<Ly;j++){
    for(i=0;i<Lx;i++){ // Loop on pixels of this object
      x=(float)(xmin+i);
      y=(float)(ymin+j);
      r2=(x-x0)*(x-x0)+(y-y0)*(y-y0);
      if ((int)r2<=(R1+1)*(R1+1)){
	s+=pixels[j*Lx+i]-bg;
	if ((int)r2>R1*R1){
	  fprintf(log1,"%f ",pixels[j*Lx+i]+1.0);
	}
	else{fprintf(log1,"%f ",pixels[j*Lx+i]);}}
      else{fprintf(log1,"%f ",pixels[j*Lx+i]);}
    }
    fprintf(log1,"\n");
  }
  printf("%d %f %f %f\n",id,R1,dR1,s);
#endif

  
  return(R1);
}
