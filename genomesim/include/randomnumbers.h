#ifndef __RANDOM_NUMBERS_H__
#define __RANDOM_NUMBERS_H__

#include <stdio.h>
#include <stdlib.h>

#include "globals.h"

#ifdef DEBUG
#define DEBUG_PRINT(...) do{ fprintf(stderr,__VA_ARGS__);}while(0)
#else
#define DEBUG_PRINT(...) do{} while(0)
#endif

#define MAX(A,B)	((A) >(B)) ? (A) : (B)

#define RND_NUM_BUF_SZ	2000

extern numeric_t_rp logOneMinusU[RND_NUM_BUF_SZ];
extern const numeric_t_rp * logOneMinusUEnd;
extern numeric_t_rp * currLogOneMinusUPos;



void randomnumbers_init(int seed);
void randomnumbers_free();
void rand_int(size_t lim,int * val);

void randomUniqueIntArray(size_t nElems,size_t nMin,size_t nMax,unsigned *dest);
void randomNonUniqueIntArray(size_t nElems,size_t nMin, size_t nMax, unsigned *dest);
void randomBinomialArray(size_t len, size_t d_low,size_t d_high,numeric_t_rp p, unsigned * dest);

void randomBoolArray(size_t nElems,unsigned int nOnes,unsigned char * dest);

void randomUniformNumberArray(size_t nElems,numeric_t_rp min, numeric_t_rp max,numeric_t_rp *dest);
rate_t_rp getNextExponential(rate_t_rp rate);
unsigned getNextBindingPosition(unsigned max);

// size_t sumArray(size_t len, int *arr);


void shufflePtrArray(void * data, const size_t len, const size_t sz);

void * getStream();



#endif