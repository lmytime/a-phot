#include <fitsio.h>
#define FLUXPRIOR

typedef unsigned char uchar;

typedef struct structcatalog {
  int num;
  int id;
  double x_image;
  double y_image;
  int xmin_image;
  int xmax_image;
  int ymin_image;
  int ymax_image;
  double detbckg;
  double measbckg;
#ifdef FLUXPRIOR
  double fpf,pf,epf;
#endif
} catalogstruct;

typedef struct structfilter {
  double *conv;	/* pointer to the convolution mask */
  int nconv;	/* total number of elements */
  int convw;	/* x size of the mask */
  int convh;	/* y size of the mask */
  float	varnorm;
} filterstruct;

typedef struct structthumbnail {
  int id;
  double x_image;
  double y_image;
  int xmin_image;
  int xmax_image;
  int ymin_image;
  int ymax_image;
  int width;
  int height;
  uchar position; // = 1 inside; = 0 bound
  uchar clipped; // = 1 yes; = 0 no
  unsigned int isoarea;
  unsigned int corrisoarea;
  double detnorm;
  double detbckg;
  double detbckgsigma;
  double measnorm;
  double measbckg;
  double measbckgsigma;
  double meanrms;
  double sigmarms;
  double chi2;
  double corrchi2;
  double correction;
  double fitcorrection;
  double *segmentation;
  double *fitsegmentation;
  double *model;
  double *selectedmodel;
  double *measure;
  double *rms;
  fitsfile *thsegmentation;
  fitsfile *thfitsegmentation;
  fitsfile *thmodel;
  fitsfile *thselectedmodel;
  fitsfile *thmeasure;
  fitsfile *thrms;
} thumbnailstruct;

#define OUTPARAM_BUFSIZE 10

typedef struct structoutparam {
  int id;
  double x_image;
  double y_image;
  int xmin_image;
  int xmax_image;
  int ymin_image;
  int ymax_image;
  int width;
  int height;
  unsigned int isoarea;
  unsigned int corrisoarea;
  int contaminant;
  double detnorm;
  double detbckg;
  double detbckgsigma;
  double measnorm;
  double measbckg;
  double measbckgsigma;
  double meanrms;
  char position[OUTPARAM_BUFSIZE];
  char status[OUTPARAM_BUFSIZE];
  double fflux;
  double ffluxerr;
  double detmag;
  double measmag;
  double measmagerr;
  double modelmag;
  double modelmagerr;
  double fluxratio;
  double chi2;
  double redchi2;
  double corrchi2;
  double redcorrchi2;
  double correction;
  double fitcorrection;
  double b;
  double bfit;
  double bdiff;
} outparamstruct;

typedef struct structbckgparam {
  int id;
  unsigned long numpix;
  unsigned long niter;
  double mean;
  double median;
  double mode;
  double sigma0;
  double sigmavar;
  double bckg;
  double bckgsigma;
  char type[OUTPARAM_BUFSIZE];
} bckgparamstruct;

typedef struct structbckgoutparam {
  int id;
  int detdlt;
  int detndlt;
  unsigned long detnumpix;
  unsigned long detniter;
  double detmean;
  double detmedian;
  double detmode;
  double detsigma0;
  double detsigmavar;
  double detbckg;
  double detbckgsigma;
  char dettype[OUTPARAM_BUFSIZE];
  int measdlt;
  int measndlt;
  unsigned long measnumpix;
  unsigned long measniter;
  double measmean;
  double measmedian;
  double measmode;
  double meassigma0;
  double meassigmavar;
  double measbckg;
  double measbckgsigma;
  char meastype[OUTPARAM_BUFSIZE];
} bckgoutparamstruct;

typedef struct structalignoutparam {
  int id;
  double x_image;
  double y_image;
  double detnorm;
  double measnorm;
  double detx;
  double dety;
  double measx;
  double measy;
  unsigned long niter;
  double offx;
  double offy;
  double shiftx;
  double shifty;
  long roundoffx;
  long roundoffy;
  double detSN;
  double measSN;
  char alignment[OUTPARAM_BUFSIZE];  
} alignoutparamstruct;
