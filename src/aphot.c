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


/*
  A-PHOT: aperture photometry package
  Computes aperture photometry on chosen image image.
  E. Merlin 2015-...
*/

//#define VERBOSE
//#define DEBUG
#define SUBGRID
#define WRITEAPERIM
//#define POLYPA
#define MAG_CAT


double *pixels, *pixels_orig, *pixels_rms, *pix_aper, *fff, **dw;
int *pixels_seg, *pixels_seg_orig, *pixels_flag;
long fpixel[2], lpixel[2]; //, inc[2] = {1,1};
long naxes[2] = {1,1};
fitsfile *fptr, *fptrrms, *fptrseg, *fptrflag;
float CRPIX1ref=0.0,CRPIX2ref=0.0;
float CRPIX1=0.0,CRPIX2=0.0;
double CRVAL1=0.e0, CRVAL2=0.e0;
double CD1_1=0.0; double CD2_2=0.0;
char *CTYPE1, *CTYPE2;
double TOL = 1.e-12;
double RMSTOL = 1000.;
long MeasWidth, MeasHeight;
FILE *olog,*ulog,*mlog,*llog,*CatalogFile,*AperFile;
char buffer[BUFFERSIZE];
double dx, dy, dsubpix, dsx, dsy;
int iii, jjj, i, j, W, H;
int find_best_sn=0, max_iter_sn=6;
int bkgd=0, symreplace=0, obj_flag, pix_max, blend_obj;
float a,b,theta;
float bad_pix, contamin_pix, aper_area;
int xc, yc, xmin, ymin, xmax, ymax;
int Area, Lx, Ly;
float fac_kron=1.0;

float subgrid(double*, double **, int, int, float, float, float, float, int, double, double*, int, int, double, float, float, double, double);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
void writeaper(double*, int, int, char*);
double kron(int, int, int, int, float, float, float, float, float, double*, int*, int*, float, int, float, float);
double background(double*, double*, int*, int*, int, int, int, float, float, int, int, int, int, int, float, float, int, int);
extern void readpixels(int, int, int, int, int);
extern void freepixels();
extern void moments(int, int, int, int, float, float, double*, double*, int*, int*, float, int, float);
extern void replace_sym(int, int, int, int, int, double);
extern int phot_flag(float, float, float, double, int, int, float, float, double, double, double, double, int, int);


int main(int argc, char *argv[])
{
  int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
  int bitpix, naxis;
  int id, rrr, rbuf, px, py, pxmin, pxmax, pymin, pymax;
  int NumProcessedSources=0;
  float x0, y0, x1, y1, x, y, xp, yp, xm, ym, xa, ya, xb, yb;
  float e, aaa, bbb, aaa0, bbb0, pa; 
  double PA, cosPA, sinPA;
  int naper, keep_pix, naperk, ni, nj;
  double *aper, *aperk, facarea;
  double R, R2, Rmax, RRR, fmax;
  double *f, *errf, *fell, *errfell, fseg, errfseg;
  double fbkgd, mumax;
  char *logbasename, *logbasenamemorph; 
  char logname[BUFFERSIZE];
  int binsubpix, ap, apk, useseg, useflag, doclip, areaseg;
  float num, GAIN;
  float fac_rbkgd, fac_sigma, fac_sigma2, fmin_kron, fmax_kron;
  float v_aper=1000.0;
  float R1, ZP=0.0, pixsc=1.0, mjyfac=0.0;
  int clipmax=0;

  // Parameters
  fac_rbkgd=atof(argv[18]);
  rbuf=atoi(argv[17]);
  fac_sigma2=atof(argv[16]);
  fac_sigma=atof(argv[15]);
  fmax_kron=atof(argv[14]);
  fmin_kron=atof(argv[13]);
  printf("Parameters: fac_rbkgd, rbuf, fac_sigma, fac_sigma2, fmin_kron, fmax_kron\n");
  printf("%f %d %f %f %f %f\n",fac_rbkgd, rbuf, fac_sigma, fac_sigma2, fmin_kron, fmax_kron);

  binsubpix=atoi(argv[7]);
  dsubpix=1.0/(float)binsubpix;
  
  ZP=atof(argv[19]);
  RMSTOL=atof(argv[23]);
  
  mjyfac=atof(argv[24]);
  if (mjyfac>0.0)
    {mjyfac=pow(10,-0.4*(ZP-23.9));
      printf("Saving catalog in microJy\n");}
  else
    {mjyfac=1.0;}

  int force_kron_comput=0;
  force_kron_comput=atoi(argv[25]);
  if (force_kron_comput>0) printf("Forcing computation of Kron radius\n");
  
  // Compute subgrid distances
  fff = (double *) malloc(binsubpix*binsubpix * sizeof(double));
  dw = dmatrix(0,8,0,(long)(binsubpix*binsubpix-1));

  for (jjj=0;jjj<binsubpix;jjj++)
    for (iii=0;iii<binsubpix;iii++)
      {
	for (nj=0;nj<3;nj++)
	  for (ni=0;ni<3;ni++)
	    {
	      // compute weight as the inverse of distance from the center of this pixel
	      dw[nj*3+ni][jjj*binsubpix+iii]=1./MAX(0.01,sqrt(pow(((1.0+((float)iii+0.5)*dsubpix)-(ni+0.5)),2)+pow(((1.0+((float)jjj+0.5)*dsubpix)-(nj+0.5)),2)));
	    }
      }

  // Open image
  printf("Reading Image %s\n",argv[1]);
  if(fits_open_file(&fptr, argv[1], READONLY, &status)){
    fits_report_error(stderr, status);
    printf("ERROR: No image found. Aborting\n");
    return(RETURN_FAILURE);
  }
  else{
    /* Read size */
    if(fits_get_img_size(fptr, 2, naxes, &status)) {
      fits_report_error(stderr, status);
      return(RETURN_FAILURE);  
    }
    MeasWidth=naxes[0];
    MeasHeight=naxes[1];
  }   

  W=MeasWidth;
  H=MeasHeight;

#ifdef WRITEAPERIM
  pix_aper = (double *) malloc(W*H * sizeof(double));
  for (i=0;i<W*H;i++) pix_aper[i]=0;
#endif

  GAIN=atof(argv[8]);
  logbasename=argv[11];
  logbasenamemorph=argv[12];
  
  // Open RMS Image
  printf("Reading RMS Image %s\n",argv[2]);
  if(fits_open_file(&fptrrms, argv[2], READONLY, &status)){
    fits_report_error(stderr, status);
    printf("ERROR: No RMS image found. Aborting\n");
    return(RETURN_FAILURE);
  }
  else{
    /* Read size */
    if(fits_get_img_size(fptrrms, 2, naxes, &status)) {
      fits_report_error(stderr, status);
      return(RETURN_FAILURE);  
    }
    if (W!=naxes[0] || H!=naxes[1]) 
      {
	printf("ERROR: RMS has width different from measure image\n");
	fits_report_error(stderr, status);
	return(RETURN_FAILURE);
      }
  }  

  // Open Segmentation Image
  printf("Reading Segmentation Image %s\n",argv[3]);
  useseg=1;
  if(fits_open_file(&fptrseg, argv[3], READONLY, &status)){
    printf("WARNING: No segmentation image found (status %d)\n",status);
    printf("Will proceed without segmentation info\n");
    useseg=0;
    status=0;
  }
  else{
    if(fits_get_img_size(fptrseg, 2, naxes, &status)) {
      fits_report_error(stderr, status);
      return(RETURN_FAILURE);  
    }
    if (W!=naxes[0] || H!=naxes[1]) 
      {
	printf("ERROR: SEG has width different from measure image\n");
	fits_report_error(stderr, status);
	return(RETURN_FAILURE);	
      }
  }

  // Open Flag Image
  printf("Reading Flag Image %s\n",argv[21]);
  useflag=1;
  if(fits_open_file(&fptrflag, argv[21], READONLY, &status)){
    printf("WARNING: No flag image found (status %d)\n",status);
    printf("Will proceed without flag info\n");
    useflag=0;
    status=0;
  }
  else{
    if(fits_get_img_size(fptrflag, 2, naxes, &status)) {
      fits_report_error(stderr, status);
      return(RETURN_FAILURE);  
    }
    if (W!=naxes[0] || H!=naxes[1]) 
      {
	printf("ERROR: FLAG has width different from measure image\n");
	fits_report_error(stderr, status);
	return(RETURN_FAILURE);	
      }
  }
  
  // Read circular apertures. The file must have apertures (DIAMETERS) listed from smaller to larger
  printf("Reading apertures\n");
  if ((AperFile = fopen(argv[5], "r")) == NULL) {
    printf("Error: can't open file \"%s\" in mode 'r'\n", 
	   argv[5]);
    return(RETURN_FAILURE);
  }
  ap=0; 
  while(fscanf(AperFile, "%f", &num) > 0.0) {
    ap++;
  }
  naper=ap;
  printf("Computing %d circular apertures\n",naper);
  aper = (double *) malloc(naper * sizeof(double));
  f = (double *) malloc(naper * sizeof(double)); 
  errf = (double *) malloc(naper * sizeof(double));
  rewind(AperFile);
  ap=0;
  while(fscanf(AperFile, "%f", &num) > 0.0) {
    aper[ap] = 0.5*num; // HERE change diameters to radii
    ap++;
  }
  fclose(AperFile);
  Rmax=aper[naper-1]; // Rmax is largest circular aperture for now.
  
  // Read elliptical factors.
  // The file must have factors listed from smaller to larger
  printf("Reading Kron multiplicative factors\n");
  if ((AperFile = fopen(argv[6], "r")) == NULL) {
    printf("Error: can't open file \"%s\" in mode 'r'\n", 
	   argv[6]);
    return(RETURN_FAILURE);
  }
  apk=0; 
  while(fscanf(AperFile, "%f", &num) > 0.0) {
    apk++;
  }
  naperk=apk;
  aperk = (double *) malloc(naperk * sizeof(double));
  for (i=0;i<naperk;i++) aperk[i]=0.0;
  rewind(AperFile);
  apk=0;
  while(fscanf(AperFile, "%f", &num) > 0.0) {
    aperk[apk] = num;
    apk++;
  }
  fclose(AperFile);
  
  if (naperk==1 && fabs(aperk[0])<0.001){
    printf("Searching for best S/N elliptical aperture\n");
    find_best_sn=max_iter_sn; // put counter equal to maximum (>0) to flag the option, will be decreased in the loop
    free(aperk);
    naperk=5; // 5 apertures will be needed
    aperk = (double *) malloc(naperk * sizeof(double));
    fmax=2.2; // maximum allowed value
  }
  else{
    printf("Computing %d elliptical apertures\n",naperk);
    fmax=aperk[naperk-1];
  }
  // fmax is largest factor for elliptical apertures (multiplies the value in input catalog)
  
  fell = (double *) malloc(naperk * sizeof(double)); 
  errfell = (double *) malloc(naperk * sizeof(double));
  
  // Open output files
  snprintf(logname,BUFFERSIZE,"%s",logbasename);
  if ((olog = fopen(logname, "w")) == NULL) {
    fprintf(stderr, "Error: can't open file \"%s\" in mode 'w'\n", 
	    logname);
    return(RETURN_FAILURE);
  }  
  fprintf(olog, "#id x y area_seg f_seg ");
  for (ap=0;ap<naper;ap++)
    {
      fprintf(olog, "ApC%d ",ap+1);
    }
  
  if (find_best_sn>0)
    {
      ulog = fopen("best_SN_aper.log","w");
      fprintf(ulog, "#ID fac_bestSN a_bestSN[pix] bestSN\n");
      fprintf(olog, "ApEll_bestSN ");
    }
  else
    {
      for (apk=0;apk<naperk;apk++)
	{
	  fprintf(olog, "ApEll%d ",apk+1);
	}
    }

  fprintf(olog, "errf_seg ");
  for (ap=0;ap<naper;ap++)
    {
      fprintf(olog, "err_ApC%d ",ap+1);
    }
  
  if (find_best_sn>0)
    {
      fprintf(olog, "err_ApEll_bestSN ");
    }
  else
    {
      for (apk=0;apk<naperk;apk++)
	{
	  fprintf(olog, "err_ApEll%d ",apk+1);
	}
    }

  fprintf(olog, "mu_max phot_flag ");
  
  snprintf(logname,BUFFERSIZE,"%s",logbasenamemorph);
  mlog = fopen(logname, "w");
  fprintf(mlog,"#id a e theta R1\n");

#ifdef MAG_CAT
  pixsc=atof(argv[20]);
  if (ZP>0.0){
    printf("Computing magnitudes using ZP=%f and pixelscale=%f\n",ZP,pixsc);
    snprintf(logname,BUFFERSIZE,"%s_mag",logbasename);
    llog = fopen(logname, "w");
    fprintf(llog,"#id x y area_seg mag_seg ");
    for (ap=0;ap<naper;ap++)
      {
	fprintf(llog, "ApC%d ",ap+1);
      }
    for (apk=0;apk<naperk;apk++)
      {
	fprintf(llog, "ApEll%d ",apk+1);
      }
    fprintf(llog,"errmag_seg ");
    for (ap=0;ap<naper;ap++)
      {
	fprintf(llog, "err_ApC%d ",ap+1);
      }
    for (apk=0;apk<naperk;apk++)
      {
	fprintf(llog, "err_ApEll%d ",apk+1);
      }
    fprintf(llog, "mu_max phot_flag\n");
  }
#endif
  
  bkgd=atoi(argv[9]);
  
  if (bkgd>0) {
    printf("Measuring and subtracting local background\n");
    clipmax=20;
    if (useseg==0) clipmax=10;
    fprintf(olog, "bkgd_pix ");
    //fmax=MAX(fmax,2.0); // If needed enlarge to size needed to measure background
  }
  fprintf(olog, "\n");
  
  doclip=atoi(argv[10]);
  if (doclip>0) {
    printf("Clipping out pixels below -1.0sigma\n"); //,fac_sigma2);
  }

  symreplace=atoi(argv[22]);
  if (symreplace>0){
    printf("Replacing bad/contaminated pixels with values from symmetric ones\n");
  }
  
  // Read in catalog and do everything :)
  printf("Reading catalog\n");
  if ((CatalogFile = fopen(argv[4], "r")) == NULL) {
    fprintf(stderr, "Error: can't open file \"%s\" in mode 'r'\n", 
	    argv[4]);
    return(RETURN_FAILURE);
  }
  
  printf("Working...\n"); 
  
  while (fgets(buffer, BUFFERSIZE, CatalogFile) != NULL && !feof(CatalogFile))  // reads the catalog
    { // Loop on objects in input catalog
      
      // Skip line if it begins with "#" or the input columns are not 6 or ID is < 0 
      // NOTE: a is the major semiaxis of the aperture, not of the object!
      if ((buffer - '#') == 0 || 
	  (sscanf(buffer, "%d\t%g\t%g\t%g\t%g\t%g", &id, &x0, &y0, &a, &e, &pa)) != 6 ||
	  (id < 1)) {
	continue;}    
      
#ifdef VERBOSE
      printf(">>> Read in values for ID %d: %f %f %f %g %g\n",id,a,e,b,pa,PA);
#endif
      
      NumProcessedSources++;
      
#ifdef VERBOSE
      printf(">Working on ID %d\n",id);
#else
      printf("\33[1M>Working on ID %d\n\33[1A",id);     
#endif
      
      // C arrays start from zero
      x0-=1.0;
      y0-=1.0;
      
      // central pixel
      xc=floor(x0+0.5e0);
      yc=floor(y0+0.5e0);
      
      // Now must determine the region of interest. It has a radius given by
      // rrr = MAX (fmax*R1, R1+rbuf, Rmax, r_segm) where:
      // R1 is the Kron radius, either entered as input or computed in the following,
      // so first number is for elliptical aperture photometry measurement and
      // second one is for background measurement (if computed);
      // Rmax is the largest circular aperture input by the user so
      // third number is for circular aperture photometry; and
      // r_segm is a dimension input by the user to internally compute
      // morphological parameters
      
      rrr=(int)ceil(Rmax); // so start with the circular apertures
      
#ifdef VERBOSE
      printf("Reading in pixels values\n");
#endif

      mumax=0.0;
      fbkgd=0.0;
      if (a>0.0)
	// Parameters of the ellipse are input by the user
	{

	  b=a*(1.0-e);
	  rrr=MAX(rrr,(int)ceil(fmax*a)); // this for elliptical apertures
	  rrr=MAX(rrr,(int)ceil(fac_rbkgd*a+rbuf)); // this is for background
	  // that's it - no need to compute parameters internally so r_segm does not matter
	  
	  readpixels(xc,yc,rrr+1,useseg,useflag);
	  if (symreplace>0) replace_sym(xc-xmin,yc-ymin,id,Lx,Ly,RMSTOL);

	  // Estimate background (or do clipping)
	  if (bkgd+doclip>0) fbkgd=background(pixels, pixels_rms, pixels_seg, pixels_flag, Lx, Ly, fac_rbkgd*a, x0, y0, xmin, ymin, clipmax, rbuf, id, fac_sigma, fac_sigma2, doclip, bkgd);

	  if (force_kron_comput==1)
	    {
	      // read in parameters of the ellipse but recompute Kron radius
	      R1=fac_kron*kron(Lx,Ly,xmin,ymin,x0,y0,a,b,aper[naper-1],pixels,pixels_seg,pixels_flag,fbkgd,id,fmin_kron,fmax_kron);
	      fprintf(mlog,"%d %f %f %f %f\n",id,a,e,theta,R1);
	      a=R1;
	      b=a*(1.0-e);
	    }
	  
	}
      else
	{
	  // Will estimate morpho parameters internally
	  // -a is r_segm, a radius of the region including the object segmentation
	  // or estimated with some XMIN-XMAX YMIN-YMAX method.

	  if (useseg==0)
	    {
	      printf("TERROR: (currently) need segmentation to compute morphological paramenters. Aborting\n");
	      return(1);
	    }
	  
	  rrr=MAX(rrr,(int)ceil(-a)); // start using the dimension input by the user 
	  
	  int read_pix=1;
	  while(read_pix>0)
	    {
#ifdef VERBOSE
	      printf("_______ iter %d, radius %d, useseg: %d, useflag: %d\n",read_pix,rrr+rbuf,useseg,useflag);
#endif
	      // Now read in pixels values to compute morpho parameters
	      if (read_pix<3)
		{readpixels(xc,yc,rrr+1,useseg,useflag);
		  if (symreplace>0) replace_sym(xc-xmin,yc-ymin,id,Lx,Ly,RMSTOL);
		}

	      // Compute morpho parameters a,b,theta
	      moments(Lx, Ly, xmin, ymin, x0, y0, pixels, pixels_rms, pixels_seg,pixels_flag,fbkgd,id,pa);
	      pa=theta*57.2958;
	      // Estimate Kron radius
	      R1=fac_kron*kron(Lx,Ly,xmin,ymin,x0,y0,a,b,aper[naper-1],pixels,pixels_seg,pixels_flag,fbkgd,id,fmin_kron,fmax_kron);
	      e=1-b/a;
	      
#ifdef VERBOSE
	      printf("Kron: id %d a %f e %f theta %f R1 %f\n",id,a,e,theta*57.2958,R1);
#endif

	      if ((int)ceil(MAX(fmax*R1,fac_rbkgd*R1+rbuf))>rrr) // the computed R1 is larger than the region: re-read in pixel values
		{
#ifdef VERBOSE
		  printf("Enlarging and reading again pixels values: r= %d -> %d\n",rrr,(int)ceil(MAX(fmax*R1,R1+rbuf)));
#endif
		  rrr=(int)ceil(MAX(fmax*R1,fac_rbkgd*R1+rbuf));
		  read_pix=1;
		  freepixels();
		}
	      
	      else
		{
		  // Estimate background (first round using the first guess, then with final morpho params)
		  if (bkgd+doclip>0)
		    {
#ifdef VERBOSE
		      printf(">> computing bkg\n");
#endif
		      fbkgd=background(pixels, pixels_rms, pixels_seg, pixels_flag, Lx, Ly, fac_rbkgd*R1, x0, y0, xmin, ymin, clipmax, rbuf, id, fac_sigma, fac_sigma2, doclip, bkgd);
		    }

		  
		  if (read_pix>1)
		    {
		      read_pix=0;
		    } 
		  else
		    {
		      if ((int)ceil(MAX(fmax*R1,R1+rbuf))>rrr)
			{
#ifdef VERBOSE
			  printf("Enlarging and reading again pixels values: r= %d -> %d\n",rrr,(int)ceil(MAX(fmax*R1,R1+rbuf)));
#endif
			  rrr=(int)ceil(MAX(fmax*R1,R1+rbuf));
			  read_pix=2;
			  freepixels();
			}
		      else
			{
			  read_pix=3;
			}
		    }
		  /*
		  if (read_pix==1)
		    {
		      read_pix=0;
		    }
		  */
		}
	    }
	  
	  fprintf(mlog,"%d %f %f %f %f\n",id,a,e,theta*57.2958,R1);
	  a=R1; // R1 is Kron Radius; set the semi-axis of the elliptical aperture equal to it
	  b=a*(1.0-e);
	  
	} 
	      
      /*
      for (j=0;j<Ly;j++){
	for (i=0;i<Lx;i++){
	  printf("%f ",pixels[j*Lx+i]);
	}
	printf("\n");
      } 
      */

      // Position angle stuff
      
      PA=(double)(pa*2.0e0*3.14159265358979323846/360.0e0); // deg to rad

#ifdef POLYPA
      // From http://lab.polygonal.de/?p=205
      //always wrap input angle to -PI..PI
      if (PA < -3.14159265358979323846)
	PA += 6.283185307179586;
      else
	if (PA >  3.14159265358979323846)
	  PA -= 6.283185307179586;
      //compute sinPAe
      if (PA < 0)
	{
	  sinPA = 1.27323954 * PA + .405284735 * PA * PA;    
	  if (sinPA < 0)
	    sinPA = .225 * (sinPA *-sinPA - sinPA) + sinPA;
	  else
	    sinPA = .225 * (sinPA * sinPA - sinPA) + sinPA;
	}
      else
	{
	  sinPA = 1.27323954 * PA - 0.405284735 * PA * PA;    
	  if (sinPA < 0)
	    sinPA = .225 * (sinPA *-sinPA - sinPA) + sinPA;
	  else
	    sinPA = .225 * (sinPA * sinPA - sinPA) + sinPA;
	}
      //compute cosine: sin(PA + PI/2) = cos(PA)
      PA += 1.5707963267948966;
      if (PA >  3.14159265358979323846)
	PA -= 6.283185307179586;
      if (PA < 0)
	{
	  cosPA = 1.27323954 * PA + 0.405284735 * PA * PA;
    
	  if (cosPA < 0)
	    cosPA = .225 * (cosPA *-cosPA - cosPA) + cosPA;
	  else
	    cosPA = .225 * (cosPA * cosPA - cosPA) + cosPA;
	}
      else
	{
	  cosPA = 1.27323954 * PA - 0.405284735 * PA * PA;
	  if (cosPA < 0)
	    cosPA = .225 * (cosPA *-cosPA - cosPA) + cosPA;
	  else
	    cosPA = .225 * (cosPA * cosPA - cosPA) + cosPA;
	}
#else
      double DPA=2.0*PA;
      double cos2PA;
      if (DPA>2.0e0*3.14159265358979323846) DPA-=2.0e0*3.14159265358979323846;
      cosPA=cosl(PA);
      sinPA=sinl(PA);
      cos2PA=cosl(DPA);
#endif     

      for(ap=naper-1;ap>-1;ap--) // Initialize circ apertures
	{
	  f[ap]=0.0 ; 
	  errf[ap]=0.0; 
	}
      if (find_best_sn>0) // Initialize elliptical apertures multiplicative factors for best SN computation
	{
	  aperk[0]=0.1;	  
	  for (apk=1;apk<naperk;apk++)
	    {
	      aperk[apk]=aperk[apk-1]+0.2; // So initial factors are 0.1,0.3,0.5,0.7,0.9
	    }
	  
	}

      obj_flag=0;
      pix_max=0;
      blend_obj=0;
      bad_pix=0.0;
      contamin_pix=0.0;
      aper_area=3.14*a*b;
      fseg=0.0;
      errfseg=0.0;
      areaseg=0;
      
      do { // Loop for best S/N bisection (will enter only once if not needed)

	for(apk=naperk-1;apk>-1;apk--) // Initialize elliptical apertures
	  {
	    fell[apk]=0.0 ; 
	    errfell[apk]=0.0; 
	  }

	for (j=0;j<Ly;j++)
	  for(i=0;i<Lx;i++) // Loop on pixels of this object
	    {
	      x=(float)(xmin+i);
	      y=(float)(ymin+j);
	      xp=x+0.5;
	      yp=y+0.5;
	      xm=x-0.5; // lower x bound of i pixel
	      ym=y-0.5; // lower y bound of i pixel

	      if (pixels_seg_orig[j*Lx+i]==0 || pixels_seg_orig[j*Lx+i]==id)
		{
		  
		  if (fabs(pixels_orig[j*Lx+i])>TOL && pixels_rms[j*Lx+i]<RMSTOL && pixels_flag[j*Lx+i]==0)
		    {
		      
		      // "a"=farthest point; "b"=closest point
		      if (y>(float)yc) // Upper region
			{
			  if (x>(float)xc){xa=xp; ya=yp; xb=xm; yb=ym;} // 1st quad
			  else if (x<(float)xc){xa=xm; ya=yp; xb=xp; yb=ym;} // 2nd quad
			  else { // Same column, above
			    if (x0>=xc) {xa=xm; ya=yp; xb=x0; yb=ym;}
			    else {xa=xp; ya=yp; xb=x0; yb=ym;}
			  }
			}
		      else if (y<(float)yc) // Lower region
			{
			  if (x<(float)xc){xa=xm; ya=ym; xb=xp; yb=yp;} // 3rd quad
			  else if (x>(float)xc){xa=xp; ya=ym; xb=xm; yb=yp;} // 4th quad
			  else{ // Same column, below
			    if (x0>=xc) {xa=xm; ya=ym; xb=x0; yb=yp;}
			    else {xa=xp; ya=ym; xb=x0; yb=yp;}
			  }
			}
		      else // Same row
			{
			  if (x>xc){ // Right
			    if (y0>=yc) {xa=xp; ya=ym; xb=xm; yb=y0;}
			    else {xa=xp; ya=yp; xb=xm; yb=y0;}
			  }
			  else if (x<xc){ // Left
			    if (y0>=yc) {xa=xm; ya=ym; xb=xp; yb=y0;}
			    else{xa=xm; ya=yp; xb=xp; yb=y0;}
			  }
			  else { // Same pixel!
			    xb=x0; yb=y0; 
			    if(x0>xc){
			      xa=xm;}
			    else{
			      xa=xp;}
			    if(y0>yc){
			      ya=ym;}
			    else{
			      ya=yp;}
			  }
			}
		      
		      if (find_best_sn==0 || find_best_sn==max_iter_sn) { // This is for circular apertures, so: only enter once, i.e. if no loop for bestSN is needed or if it's the 1st time here
			
			keep_pix=1;

			if (pixels_seg_orig[j*Lx+i]==id)
			  {
			    fseg+=(pixels[j*Lx+i]-fbkgd);
			    errfseg+=pixels_rms[j*Lx+i]*pixels_rms[j*Lx+i];
			    areaseg+=1;
			  }			
			
			for(ap=naper-1;ap>-1;ap--) // Loop on circular apertures, starting from larger
			  {
			    facarea=0.0;
			    if (keep_pix) // If the pixel has been excluded before, no need to check it again
			      {
				R=aper[ap]; // First loop has R=rrr
				R2=R*R;
				facarea=1.0;
				
				// Check if the pixel is totally outside R
				if ((xb-x0)*(xb-x0)+(yb-y0)*(yb-y0)>R2){
				  facarea=0.0;
				  // Also, don't need to check this pixel anymore
				  keep_pix=0;
				}
				// Now check if it's totally inside R
				else if ((xa-x0)*(xa-x0)+(ya-y0)*(ya-y0)>=R2){ // Nope, do subgrid			  
#ifdef SUBGRID
				  facarea=subgrid(fff,dw,i,j,xm,ym,x0,y0,binsubpix,dsubpix,pixels,Lx,Ly,R2,-1.0,0.0,0.0,0.0);
#else
				  facarea=0.0;
#endif
				}
				else // well it seems it's inside
				  {
				    if (pixels[j*Lx+i]-mumax>0.0001)
				      {
					mumax=pixels[j*Lx+i];
					pix_max=0;
				      }
				    if (fabs(pixels[j*Lx+i]-mumax)<0.0001*mumax)
				      {
					pix_max+=1;
				      }
				  }
			      }
			    
			    f[ap]+=facarea*pixels[j*Lx+i];
			    errf[ap]+=(facarea*pixels_rms[j*Lx+i])*(facarea*pixels_rms[j*Lx+i]);
			    
			    // Subtract background?
			    f[ap]-=facarea*fbkgd;
#ifdef VERBOSE
			    //printf("o-> %d %d %d %g %g %g\n",xmin+i,ymin+j,ap,f[ap],facarea,pixels[j*Lx+i]);
#endif		    
			  } // end loop on aper
			
		      } // end IF for best S/N bisection check
		      
		      
		      // Now do elliptical apertures
		      
		      keep_pix=1;
		      
		      for(apk=naperk-1;apk>-1;apk--) // Loop on apertures, starting from larger. These apertures can be defined by the user, or they can be the 5 apertures for best S/N
			{
			  
			  facarea=0.0;
			  if (keep_pix) // If the pixel has been excluded before, no need to check it again
			    {
			      aaa=aperk[apk]*a; // semi-axis of this aperture in pixels; it is factor*a, which can be given by the user or be the Kron radius computed internally
			      bbb=aperk[apk]*b;
			      
			      facarea=1.0;
			      x1=(xb-x0)*cosPA+(yb-y0)*sinPA;
			      y1=-(xb-x0)*sinPA+(yb-y0)*cosPA;
			      if (((x1*x1)/(aaa*aaa) + (y1*y1)/(bbb*bbb)) > 1.0) // pixel totally outside ellipse
				{
				  facarea=0.0;
				  keep_pix=0.0;
				}
			      else { // now check if it's totally inside
				x1=(xa-x0)*cosPA+(ya-y0)*sinPA;
				y1=-(xa-x0)*sinPA+(ya-y0)*cosPA;
				if (((x1*x1)/(aaa*aaa) + (y1*y1)/(bbb*bbb)) >= 1.0) // nope, do subgrid
				  {
#ifdef SUBGRID
				    facarea=subgrid(fff,dw,i,j,xm,ym,x0,y0,binsubpix,dsubpix,pixels,Lx,Ly,R2,aaa,bbb,cosPA,sinPA);
#else
				    facarea=0.0;
#endif
				    
#ifdef WRITEAPERIM
				    //printf("%d %d %g\n",xmin+i,ymin+j,PA);
				    pix_aper[(ymin+j)*W+(xmin+i)]=v_aper;
#endif
				  }			    
			      }
			      
			      fell[apk]+=facarea*pixels[j*Lx+i];
			      errfell[apk]+=(facarea*pixels_rms[j*Lx+i])*(facarea*pixels_rms[j*Lx+i]);
			      
#ifdef DEBUG
			      printf("--> %d %g %g %g %g\n",apk,fell[apk],errfell[apk],facarea,pixels[j*Lx+i]);
#endif		
			      // Subtract background?
			      fell[apk]-=facarea*fbkgd;		
			      
			    } // end check keep_pix
			  
			} // end loop on elliptical aper

		    }
		  else
		    {
#ifdef DEBUG
		      printf("xxx BAD PIXEL %d %d\n",i,j);
#endif
		      // this is  a "bad" pixel: either the RMS is too high or it is flagged
		      x1=(x-x0)*cosPA+(y-y0)*sinPA;
		      y1=-(x-x0)*sinPA+(y-y0)*cosPA;
		      if (((x1*x1)/(a*a) + (y1*y1)/(b*b)) <= 1.0)
			{
			  bad_pix+=1;
			}		      
		      //if (pixels_seg_orig[j*Lx+i]==id) {bad_pix+=1;}
		    }
		}
	      else
		{
		  // this pixels is assigned to another object
#ifdef DEBUG
		  printf("xxx BLENDING/CONTAMINATION %d %d\n",i,j);
#endif		  
		  
		  // first check if it is blended
		  if (blend_obj==0)
		    {
		      pxmin=-1;
		      if (i-pxmin<0) pxmin=0;
		      pymin=-1;
		      if (j-pymin<0) pymin=0;
		      pxmax=+1;
		      if (i+pxmax>=Lx) pxmax=0;
		      pymax=+1;
		      if (j+pymax>=Ly) pymax=0;		      
		      
		      for (px=-pxmin;px<=+pxmax;px++)
			for (py=-pymin;py<=+pymax;py++)
			  {
			    if (pixels_seg_orig[(j+py)*Lx+(i+px)]==id)
			      {
				blend_obj+=1;
				break;
			      }
			  }
		    }
		  
		  // then check if it contaminates aperture flux
		  x1=(x-x0)*cosPA+(y-y0)*sinPA;
		  y1=-(x-x0)*sinPA+(y-y0)*cosPA;
		  if (((x1*x1)/(a*a) + (y1*y1)/(b*b)) <= 1.0)
		    {
		      contamin_pix+=pixels_orig[j*Lx+i];
		    }
		}
	      // end check on segmentation

	    } // End loop on pixels


	// Check object flagging
	obj_flag=phot_flag(bad_pix,aper_area,contamin_pix,fell[naperk-1],blend_obj,pix_max,x0,y0,aperk[naperk-1]*a,Rmax,cosPA,sinPA,W,H);
	
      	  
	// Finalize elliptical errors
	for (apk=0;apk<naperk;apk++)
	  {
	    if (GAIN>0.0){
	      errfell[apk]=sqrt(errfell[apk]+MAX(0.0,fell[apk])/GAIN);
	    }
	    else{
	      errfell[apk]=sqrt(errfell[apk]); 
	    }
	  }

#ifdef DEBUG
	printf("ell>>>> %f %f %f %f %f - %f %f %f %f %f -- %d \n",aperk[0],aperk[1],aperk[2],aperk[3],aperk[4],fell[0]/errfell[0],fell[1]/errfell[1],fell[2]/errfell[2],fell[3]/errfell[3],fell[4]/errfell[4],find_best_sn);
#endif

	if (find_best_sn>=1) // Loop for best S/N
	  {
	    if (fell[2]/errfell[2]>MAX(fell[1]/errfell[1],fell[3]/errfell[3]))
	      {
		aperk[0]=aperk[1];
		aperk[4]=aperk[3];
	      }
	    else
	      {
		if (fell[1]/errfell[1]>fell[3]/errfell[3])
		  {
		    aperk[4]=aperk[2];
		    aperk[2]=aperk[1];	    
		  }
		else
		  {
		    aperk[0]=aperk[2];
		    aperk[2]=aperk[3];	 
		  }
	      }
	    aperk[1]=0.5*(aperk[0]+aperk[2]);
	    aperk[3]=0.5*(aperk[4]+aperk[2]);
	    //printf(">>> %d\n",find_best_sn);
	    if (find_best_sn==1) 
	      {
		find_best_sn=-1;
	      }
	    else
	      {
		find_best_sn--;
	      }
	   
	  }

      } while (find_best_sn>0); // end loop for best S/N bisection

      // Finalize circular errors
      for (ap=0;ap<naper;ap++)
	{
	  if (GAIN>0.0){
	    errf[ap]=sqrt(errf[ap]+MAX(0.0,f[ap])/GAIN);
	  }
	  else{
	    errf[ap]=sqrt(errf[ap]); 
	  }
	}
      
      // Finalize seg errors
      if (GAIN>0.0){
	errfseg=sqrt(errfseg+MAX(0.0,fseg)/GAIN);
      }
      else{
	errfseg=sqrt(errfseg); 
      }      
      // Write output for this object
      fprintf(olog, "%d %f %f %d %f ",id, x0+1.0, y0+1.0, areaseg, mjyfac*fseg);
       
#ifdef MAG_CAT
      if (ZP>0.0){
	fprintf(llog,"%d %f %f %d %f ",id, x0+1.0, y0+1.0, areaseg, -2.5*log10(fseg)+ZP);
      }
#endif      
      for (ap=0;ap<naper;ap++)
	{
	  fprintf(olog,"%f ",mjyfac*f[ap]);
#ifdef MAG_CAT
	  if (ZP>0.0){  
	    fprintf(llog,"%f ",-2.5*log10(f[ap])+ZP);
	  }
#endif
	}

      if (find_best_sn<0)
	{
	  fprintf(ulog,"%d %f %f %f\n",id,aperk[2],aperk[2]*a,fell[2]/errfell[2]);
	  fprintf(olog,"%f ",fell[2]);
#ifdef MAG_CAT
	  if (ZP>0.0){
	    fprintf(llog,"%f ",-2.5*log10(fell[2])+ZP);
	  }
#endif	  
	}
      else
	{
	  for (apk=0;apk<naperk;apk++)
	    {
	      fprintf(olog,"%f ",mjyfac*fell[apk]);
#ifdef MAG_CAT
	      if (ZP>0.0){   
		fprintf(llog,"%f ",-2.5*log10(fell[apk])+ZP);
	      }
#endif	      
	    }
	}
      
      fprintf(olog,"%f ",mjyfac*errfseg);
      
#ifdef MAG_CAT
      if (ZP>0.0){ 
	fprintf(llog,"%f ",1.0857*errfseg/fseg);
      }
#endif
      
      for (ap=0;ap<naper;ap++)
	{
	  fprintf(olog,"%f ",mjyfac*errf[ap]);
#ifdef MAG_CAT
	  if (ZP>0.0){
	    fprintf(llog,"%f ",1.0857*errf[ap]/f[ap]);
	  }
#endif	  
	}
      
      if (find_best_sn<0)
	{
	  fprintf(olog,"%f ",mjyfac*errfell[2]);
#ifdef MAG_CAT
	  if (ZP>0.0){
	    fprintf(llog,"%f ",1.0857*errfell[2]/fell[2]);
	  }
#endif		  
	  find_best_sn=max_iter_sn;
	}
      else
	{
	for (apk=0;apk<naperk;apk++)
	  {
	    fprintf(olog,"%f ",mjyfac*errfell[apk]);
#ifdef MAG_CAT
	    if (ZP>0.0){
	      fprintf(llog,"%f ",1.0857*errfell[apk]/fell[apk]);
	    }
#endif	    
	  }
	}

      fprintf(olog,"%f ",mjyfac*(mumax-fbkgd));
#ifdef MAG_CAT
      if (ZP>0.0){
	fprintf(llog,"%f ",(-2.5*log10((mumax-fbkgd)/(pixsc*pixsc))+ZP));
      }
#endif
      
      fprintf(olog,"%d ",obj_flag);
#ifdef MAG_CAT
      if (ZP>0.0){
	fprintf(llog,"%d \n",obj_flag);
      }
#endif

      if (bkgd>0) fprintf(olog,"%f ",mjyfac*fbkgd);      
      fprintf(olog,"\n");
      
      freepixels();
      
    } // End loop on objects

  printf("\33[1MClosing files...\n");

  // Close output file
  fclose(olog);
  if (find_best_sn<0) fclose(ulog);
#ifdef MAG_CAT
  if (ZP>0.0){
    fclose(llog);
  }
#endif
  fclose(mlog);
  
  free(aper); 
  free(f); 
  free(errf);
  free(aperk); 
  free(fell); 
  free(errfell);
  free(fff); 
  free_dmatrix(dw,0,8,0,binsubpix*binsubpix-1);
  if (useseg) fits_close_file(fptrseg, &status);  
  fits_close_file(fptrrms, &status); 
  fits_close_file(fptr, &status);
  
#ifdef WRITEAPERIM
  printf("Writing apertures image...\n");
  writeaper(pix_aper,W,H,argv[1]);
  free(pix_aper);
#endif
  
  printf("Done. Bye\n");

  return(status);
}

