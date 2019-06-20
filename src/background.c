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

float median(int, double*);
float mean(int, double*);

double background(double *fff, double *fffrms, int *seg, int *flag, int Lx, int Ly, int siz, float x0, float y0, int xmin, int ymin, int clipmax, int rbuf, int id, float fac_sigma, float fac_sigma2, int doclip, int bkgd)
{
  
  int TOL,i,j,r2out,r2in,clip=0,ic,jc,nnn;
  double m0,m1,s,s2,sig,f,*fff0,*fff_temp;

  ic=(int)(x0-xmin);
  jc=(int)(y0-ymin);
  TOL=MIN(20,(int)(0.25*3.14*((siz+rbuf)*(siz+rbuf)-siz*siz)));

  /*
  for (j=jc-siz;j<=jc+siz;j++){
    for (i=ic-siz;i<=ic+siz;i++){
      if ((i-ic)*(i-ic)+(j-jc)*(j-jc)<=siz*siz){
	printf("%f ",fff[j*Lx+i]);
      }
      else
	{printf("0.0 ");}
    }
    printf("\n");
  }
  */ 
    
  clip=0;
  fff_temp = (double *) malloc(Lx*Ly * sizeof(double));
  
  sig=1.e10;
  m1=0.0;
  r2out=(int)ceil((siz+rbuf)*(siz+rbuf));
  r2in=(int)ceil(siz*siz);
  
  while (clip<2)
    {
      //printf("-----------------------\n");
#ifdef VERBOSE
      printf("Bkgd iter %d - Area %d\n",clip+1,Lx*Ly);
#endif
      s=0.0;
      s2=0.0;
      nnn=0;
      for (j=0;j<Ly;j++){
	for (i=0;i<Lx;i++){
	  float wf=0.0;
	  if ((i-ic)*(i-ic)+(j-jc)*(j-jc)<=r2out){
	    if ((i-ic)*(i-ic)+(j-jc)*(j-jc)>r2in){
	      if (seg[j*Lx+i]!=id){
		if (flag[j*Lx+i]==0) {
		  f=fff[j*Lx+i];
		  if (fabs(f-m1)<=fac_sigma*fffrms[j*Lx+i]){
		    if (fabs(f-m1)<=fac_sigma2*sig){
		      fff_temp[nnn]=f;
		      s2+=f*f;
		      nnn+=1;
		      wf=f;
		    }
		  }
		}
	      }
	    }
	  }
	  //printf("%f ",wf);
	}
	//printf("\n");
      }

      if (nnn>TOL)
	{
	  fff0 = (double *) malloc(nnn * sizeof(double));
	  for (i=0;i<nnn;i++) fff0[i]=fff_temp[i];
	  
	  m1=median(nnn,fff0);      
	  s=mean(nnn,fff0);
#ifdef VERBOSE
	  printf("> px %d median %f mean %f\n",nnn,m1,s);
#endif
	  free(fff0);
	  if (clip==0)
	    {
	      m0=m1;
	      sig=sqrt(s2/(float)nnn-s*s);
	    }
	  clip+=1;
	}
      else
	{
#ifdef VERBOSE
	  printf("nnn %d (%d), fac_sigma %f -> %f\n",nnn,(int)(0.25*3.14*((siz+rbuf)*(siz+rbuf)-siz*siz)),fac_sigma,2.0*fac_sigma);
#endif
	  fac_sigma+=0.5;
	  if (fac_sigma>1.e9){
	    TOL=TOL-1;
	  }
	  else if (fac_sigma>10.0){
	    fac_sigma=1.e10;
	    fac_sigma2=1.e10;
	  }
	}
    }

  if (doclip>0)
    {
      for (j=0;j<Ly;j++)
	for (i=0;i<Lx;i++){
	  if (fff[j*Lx+i]<-MIN(fffrms[j*Lx+i],sig))
	    {
	      flag[j*Lx+i]=1;
	    }
	}
    }
  
  free(fff_temp);

  if (bkgd>0)
    {
  
      if (fabs((m1-m0)/m0)<=0.2)
	{
#ifdef VERBOSE
	  printf("--> bkgdA %f\n",s);
#endif
	  return(s);
	}
      else
	{
#ifdef VERBOSE
	  printf("--> bkgdB %f\n",2.5*m1-1.5*s);
#endif      
	  return(2.5*m1-1.5*s);
	}
    }
  else
    {return(0.0);}
  
}

//----------------------------------------------------------------------//

float median(int n, double *x) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }   

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}

//------------------------------------------------------//

float mean(int m, double *a) {
  float sum=0;
  int i;
  for(i=0; i<m; i++)
    sum+=a[i];
  return(sum/(float)m);
}

//------------------------------------------------------//
