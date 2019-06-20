
#define	QCALLOC(ptr, typ, nel)			\
  ptr = (typ *)calloc((size_t)(nel),sizeof(typ))
#define	QMALLOC(ptr, typ, nel)  \
  ptr = (typ *)malloc((size_t)(nel)*sizeof(typ))
#define	QREALLOC(ptr, typ, nel) \
  ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ))
#define	QFREE(ptr) {free(ptr); ptr = NULL;}


/*
#define	QCALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)calloc((size_t)(nel),sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}
#define	QMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QREALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)realloc(ptr, (size_t)(nel)*sizeof(typ)))) \
		   error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}

#define	QFREE(ptr) \
		{free(ptr); \
		ptr = NULL;}
*/
