#ifndef GLOBALS_RP_H__
#define GLOBALS_RP_H__
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>
#include <xmmintrin.h>
#include "precision.h"

#define MAXGENENAMELENGTH	25
#define ulong_type	int
#define ltf_rp	"%d"
#define ulong_type_max	INT32_MAX
#define qty_type    int
#define MAX_JSON_FILENAME  500

#define DROPCONST(tp,elem)   *(tp *)&elem

#define DECLAREBASICARRAYTYPESTRUCT(tp,tpname)	typedef struct tp ## Array {	\
	struct tp * data;														\
	ulong_type length;											\
} tpname ## Array;														\
typedef struct tp ## ArrayIters {									\
	struct tp * const start;												\
	struct tp * const end;													\
	struct tp * curr;														\
} tpname ## ArrayIters;													\
tpname ## Array empty ## tpname ## Array();											\
tpname ## Array tpname ## Array_alloc(ulong_type length);				\
void tpname ## Array_free(tpname ## Array * tofree);						\
void tpname ## Array_extend(tpname ## Array * toRealloc, ulong_type additionalLength); \
tpname ## ArrayIters get ## tpname ## ArrayIters(const tpname ## Array * arr);		

#define DECLAREBASICARRAYTYPETYPEDEF(tp)	typedef struct tp ## Array {	\
	tp * data;														\
	ulong_type length;											\
} tp ## Array;														\
typedef struct tp ## ArrayIters {									\
	tp * const start;												\
	tp * const end;													\
	tp * curr;														\
} tp ## ArrayIters;													\
tp ## Array empty ## tp ## Array();											\
tp ## Array tp ## Array_alloc(ulong_type length);				\
void tp ## Array_free(tp ## Array * tofree);						\
void tp ## Array_extend(tp ## Array * toRealloc, ulong_type additionalLength); \
tp ## ArrayIters get ## tp ## ArrayIters(const tp ## Array * arr);		



#define DECLAREBASICARRAYTYPEPRIMITIVE(tp,tpname)	typedef struct tpname ## Array {	\
	tp * data;														\
	ulong_type length;											\
} tpname ## Array;														\
typedef struct tpname ## ArrayIters {									\
	tp * const start;												\
	tp * const end;													\
	tp * curr;														\
} tpname ## ArrayIters;													\
tpname ## Array empty ## tpname ## Array();											\
tpname ## Array tpname ## Array_alloc(ulong_type length);				\
void tpname ## Array_free(tpname ## Array * tofree);						\
void tpname ## Array_extend(tpname ## Array * toRealloc, ulong_type additionalLength); \
tpname ## ArrayIters get ## tpname ## ArrayIters(const tpname ## Array * arr);		

#define DEFINEBASICARRAYTYPEPRIMITIVE(tp,tpname)									\
tpname ## Array empty ## tpname ## Array(){ 											\
	return (tpname ## Array){.data = NULL,.length=0};}					\
tpname ## Array tpname ## Array_alloc(ulong_type length) {				\
	tpname ## Array toRet;												\
	toRet.data = my_malloc(length * sizeof(*toRet.data));			\
	memset(toRet.data,0,length*sizeof(*toRet.data));				\
	toRet.length = length; return toRet;}							\
	tpname ## ArrayIters get ## tpname ## ArrayIters(const tpname ## Array * arr){	\
	if (arr->length) return (tpname ## ArrayIters) {.start = arr->data, 					\
	.end = arr->data + arr->length, .curr = arr->data};				\
	else															\
		return (tpname ## ArrayIters){.start = NULL,.end=NULL,.curr=NULL};\
	}													\
void tpname ## Array_free(tpname ## Array * tofree) {						\
	tpname ## ArrayIters it = get ## tpname ## ArrayIters(tofree);			\
	if(it.start){														\
		for(it.curr=it.start;it.curr!=it.end;it.curr++)					\
	 tpname ## _free(it.curr);											\
	my_free(tofree->data);											\
	tofree->data = NULL;}											\
	tofree->length = 0;	}											\
void tpname ## Array_extend(tpname ## Array * toRealloc, ulong_type additionalLength){	\
	tp * newData = NULL;											\
	ulong_type szlength = toRealloc->length * sizeof(*(toRealloc->data));\
	szlength = szlength + additionalLength * sizeof(*(toRealloc->data));		\
	newData = realloc(toRealloc->data, szlength);						\
	if(!newData) {fprintf(stderr,"FAILURE TO REALLOC %s:%d\n",__FILE__,__LINE__);	\
	exit(-5);}														\
	toRealloc->data = newData;										\
	toRealloc->length = toRealloc->length + additionalLength;		\
}	


#define DEFINEBASICARRAYTYPESTRUCT(tp)									\
tp ## Array empty ## tp ## Array(){ 											\
	return (tp ## Array){.data = NULL,.length=0};}					\
tp ## Array tp ## Array_alloc(ulong_type length) {				\
	tp ## Array toRet;												\
	toRet.data = mallocDBG(length * sizeof(*toRet.data));			\
	memset(toRet.data,0,length*sizeof(*toRet.data));				\
	toRet.length = length; return toRet;}							\
	tp ## ArrayIters get ## tp ## ArrayIters(const tp ## Array * arr){	\
	tp ## ArrayIters toRet = {.start = arr->data, 					\
	.end = arr->data + arr->length, .curr = arr->data};				\
	return toRet;}													\
void tp ## Array_free(tp ## Array * tofree) {						\
	if(tofree->data == NULL) return;								\
	tp ## ArrayIters it = get ## tp ## ArrayIters(tofree);			\
	for(it.curr=it.start;it.curr!=it.end;it.curr++)					\
	 tp ## _free(it.curr);											\
	freeDBG(tofree->data);											\
	tofree->length = 0;		}										\
void tp ## Array_extend(tp ## Array * toRealloc, ulong_type additionalLength){	\
	tp * newData = NULL;											\
	ulong_type szlength = toRealloc->length * sizeof(*(toRealloc->data));\
	szlength = szlength + additionalLength * sizeof(*(toRealloc->data));		\
	newData = realloc(toRealloc->data, szlength);						\
	if(!newData) {fprintf(stderr,"FAILURE TO REALLOC %s:%d\n",__FILE__,__LINE__);	\
	exit(-5);}														\
	toRealloc->data = newData;										\
	toRealloc->length = toRealloc->length + additionalLength;		\
}																	\

#define ARRAY_TYPE_FOREACH(iter)	for(iter.curr=iter.start;iter.curr!=iter.end;iter.curr++) \
	
#define TIMEIT(func) do{clock_t start = clock(),diff;				\
func;																\
diff = clock()-start;												\
int msec = diff * 1000 / CLOCKS_PER_SEC;							\
printf("%s: %d ms\n",#func,msec);}while(0);							\


//ASSERTION DEFINITION

#ifdef WITHASSERTIONS
#include <assert.h>
#else
#undef assert
#define assert(x) do{}while(0)
#define ASSERTIF(x) do{} while(0)
#endif

//////////////////////



////////////////////

#ifndef TXC_K_ONE_HALF
#define TXC_K_ONE_HALF  6
#endif
#ifndef TXC_PROD_EXPONENT
#define TXC_PROD_EXPONENT  1
#endif

#ifndef TXL_K_ONE_HALF
#define TXL_K_ONE_HALF  5
#endif
#ifndef TXL_PROD_EXPONENT
#define TXL_PROD_EXPONENT  1
#endif

#ifndef MESS_DECAY_K_ONE_HALF
#define MESS_DECAY_K_ONE_HALF  5
#endif
#ifndef MESS_DECAY_EXPONENT
#define MESS_DECAY_EXPONENT  1
#endif


enum ConnDegreeTp {INDEGREE = 0, OUTDEGREE = 1};

#define MAX_NCONNECTIONS	255
typedef struct occ_bits_rp_t {uint64_t x[4];} occ_bits_rp_t;

#define CODING_QTY	2   //Note these are forced
#define NONCODING_QTY	4
#define MAX_MRNA_QTY	1000




#define INITLOGGINGINTERVAL	100
#define DEFAULTLOGGINGINTERVAL	1000
#define UNIFORMSAMPLINGINTERVAL	1.0


extern time_t_rp clockTime;
extern ulong_type mallocDBGCt;
extern ulong_type freeDBGCt;


extern int networkGenSeed;
extern int simulationSeed;
extern int miRNAs_active;
extern ulong_type n_cycles;



#ifdef DEBUG
#define my_malloc(x)	mallocDB(x,__func__)
#define my_free(x)	freeDB(x,__func__)
void * mallocDB(size_t sz,const char * caller);
void freeDB(void * toFree,const char * caller);
#else
#define my_malloc(x)    _mm_malloc(x,64)
#define my_free(x)      _mm_free(x)
#endif

typedef enum species_enum {UNDEFINED=0,NOSPECIES=0,CODING=1,NONCODING=2,DNA=3,MESSENGER=4,PRODUCER_SPECIES=7,
	MICRO=8,PROTEIN=16,MODULATOR_SPECIES = 24,MESSMIR=32,TFDNA=64,BOUND_SPECIES = 96} species_t;
#ifndef DEBUG
inline char const * species_names(species_t sp) {
	switch(sp) {
		case UNDEFINED: return "UNDEFINED";
		case CODING: return "CODING";
		case DNA: return "DNA";
		case NONCODING: return "NONCODING";
		case MESSENGER: return "MESSENGER";
		case MICRO: return "MICRO";
		case PROTEIN: return "PROTEIN";
		case MESSMIR: return "MESSMIR";
		case TFDNA: return "TFDNA";
		default: fprintf(stdout,"BADSPECIESNAME");exit(0);
	}
}
#else
char * species_names(species_t sp);
#endif


typedef enum SkipDisabled {INCLUDE_DISABLED=0,SKIP_DISABLED=1} SkipDisabled;

void resetTime();



DECLAREBASICARRAYTYPEPRIMITIVE(unsigned int, UnsignedInt)
DECLAREBASICARRAYTYPEPRIMITIVE(unsigned char, UnsignedChar)
void UnsignedIntArray_sort(UnsignedIntArray arr);
UnsignedIntArray UnsignedIntArray_fromArray(unsigned * start, unsigned *end);
UnsignedIntArray UnsignedArray_cat(UnsignedIntArray arr1, UnsignedIntArray arr2);
unsigned UnsignedIntArray_total(UnsignedIntArray arr);
struct cJSON;

void UnsignedCharArray_clear(UnsignedCharArray arr);
unsigned UnsignedCharArray_total(UnsignedCharArray arr);


char * readFileToString(const char * filename, size_t * strSize);

extern long commonBirthDeaths;
extern long commonBirthDeathsLast;

typedef struct synthTotals {
	ulong_type totalCodingTranscription;
	ulong_type totalNonCodingTranscription;
	ulong_type totalTranslations;
	ulong_type totalMessDecays;
	ulong_type totalForcedMirDecays;
	ulong_type totalForcedMirUnbindings;
	ulong_type totalModDecays;
	ulong_type totalMessMirBindings;
	ulong_type totalTFDNABindings;
	ulong_type totalMessMirUnBindings;
	ulong_type totalTFDNAUnBindings;
    	ulong_type nMessTranscription;
	ulong_type nMicroTranscription;
	ulong_type nTranslation;
} synthTotals;

extern synthTotals totalSynths;

#endif
