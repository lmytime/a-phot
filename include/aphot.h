#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <types.h>
#include <alloc.h>
#include <nrutil.h>
#include <fitsio.h>

#ifndef VERSION
#define VERSION "APHOT 1.0.0"
#endif
#define SYSTEM_TYPE "pc-linux"
#define MACHINE_TYPE "i386"

#define BUFFERSIZE 2048
#define BIG_BUFFERSIZE 16 * BUFFERSIZE
#define TRUE 1
#define FALSE 0
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#define RETURN_SUCCESS 0
#define RETURN_FAILURE 1
#define CHAR_MAX_VALUE TYPE_MAXIMUM (char)
#define DSTDIR_MAX_CHAR 256

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))

/* These macros work even on ones'-complement hosts (!).
   The extra casts work around common compiler bugs.  */
#define TYPE_SIGNED(t) (! ((t) 0 < (t) -1))
#define TYPE_MINIMUM(t) (TYPE_SIGNED (t) \
			 ? ~ (t) 0 << (sizeof (t) * CHAR_BIT - 1) \
			 : (t) 0)
#define TYPE_MAXIMUM(t) (~ (t) 0 - TYPE_MINIMUM (t))



