//globals.h - General macros, constants, functions and global vars
/*
    
    Copyright (C) 2018, Russell Posner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
#ifndef GLOBALS_RP_H__
#define GLOBALS_RP_H__
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdint.h>


#define MAXGENENAMELENGTH	25
#define ulong_type	unsigned long
#define ltf_rp	"%zu"
#define ulong_type_max	UINT64_MAX
static const unsigned MAX_JSON_FILENAME = 500;

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
	toRet.data = mallocDBG(length * sizeof(*toRet.data));			\
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
	freeDBG(tofree->data);											\
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

//TESTING BIMODAL CONNECTIONS THIS IS REALLY HACKY
#ifndef BIMODAL_MIRNA_RP
#define BIMODAL_MIRNA_RP	0
#define BIMODAL_P2_RP	0
#define BIMODAL_SELECT_RP	1.0
#else
#define BIMODAL_P2_RP	0.1
#define BIMODAL_SELECT_RP	0.2
#endif

////////////////////


enum ConnDegreeTp {INDEGREE = 0, OUTDEGREE = 1};
//This is important - there are hard maxes here 
//Max # of modulators is 63 (using 64 bit units to represent occupancy), can be changed by
//changing type of bitsandmod in basicgeneticelement_p.h

#define MAX_NCONNECTIONS	63
typedef uint64_t occ_bits_rp;
//default quantities of DNA when not otherwise specified
#define CODING_QTY	2
#define NONCODING_QTY	8
//Also important - specifies max # of mRNA molecules per gene
//Can make this bigger or smaller freely - system will check to make sure this number isn't
//exceeded and will error
#define MAX_MRNA_QTY	1000

//Precision definitions
#ifdef DOUBLE_PRECISION_RP
typedef double numeric_t_rp;
typedef double rate_t_rp;
typedef double time_t_rp;
typedef double effect_t_rp;
#define MEPS_RP	DBL_EPSILON
#else
typedef float numeric_t_rp;
typedef float rate_t_rp;
typedef float time_t_rp;
typedef float effect_t_rp;
#define MEPS_RP	FLT_EPSILON
#endif

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


extern ulong_type totalProductions;
extern ulong_type totalDecays;
extern ulong_type totalBindings;
extern ulong_type totalUnbindings;
extern ulong_type totalForcedUnbindings;
extern ulong_type totalForcedDecays;

//debug versions of mem funcs
#define mallocDBG(x)	mallocDB(x,__func__)
#define freeDBG(x)	freeDB(x,__func__)
void * mallocDB(size_t sz,const char * caller);
void freeDB(void * toFree,const char * caller);

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
	ulong_type nMessTranscription;
	ulong_type nMicroTranscription;
	ulong_type nTranslation;
} synthTotals;

extern synthTotals totalSynths;

#endif
