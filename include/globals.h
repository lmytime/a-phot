
#define BIG_NUMBER 1.0e+32 /* A big number */
#define VERY_BIG_NUMBER 1.0e+64 /* A very big number */
#define HUGE_NUMBER 1.0e+128 /* A huge number */
#define SMALL_NUMBER 1.0e-32 /* A small number */
#define VERY_SMALL_NUMBER 1.0e-64 /* A very small number */
#define INFINITESIMAL_NUMBER 1.0e-128 /* A infinitesimal number */
#define INIT(name); {progname = name;}
#define I_SIGN(d) (d >= 0 ? 1 : -1)
#define HALFPI 1.570796327

catalogstruct *Catalog;
filterstruct *Filter;
thumbnailstruct *Thumbnail;
outparamstruct *OutParam;
bckgparamstruct *BckgParam;
bckgoutparamstruct *BckgOutParam;
alignoutparamstruct *AlignOutParam;
uchar interactive;
uchar subdetbckg;
uchar submeasbckg;
uchar writefits;
uchar align;
uchar minimization;
uchar loadth;
uchar resume;
uchar saveth;
uchar verbose;
int NumSources;
int NumInSources;
int NDumpedSources;
int NDumpedProcessedSources;
char *progname;
