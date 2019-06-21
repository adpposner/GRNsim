//randomnumbers.h - interface for random number & shuffle functions
//Provides access to basic MKL VSL functions
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
#ifndef __RANDOM_NUMBERS_H__
#define __RANDOM_NUMBERS_H__

#include <stdio.h>
#include <stdlib.h>
#include <mkl_types.h>
#include <mkl.h>
#include <mkl_vsl.h>
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
void randomUniformDoubleArray(size_t nElems,double min,double max,double * dest);
void randomUniformSingleArray(size_t nElems,float min,float max, float * dest);
void randomBinomialArray(size_t len, size_t d_low,size_t d_high,double p, unsigned * dest);
void randomBimodalBinomialArray(size_t len, size_t d_low, size_t d_high, double p1, double p2, double p_select,
 unsigned * dest);
void randomBoolArray(size_t nElems,unsigned int nOnes,unsigned char * dest);

void randomUniformNumberArray(size_t nElems,effect_t_rp min, effect_t_rp max,effect_t_rp *dest);
rate_t_rp getNextExponential(rate_t_rp rate);
unsigned getNextBindingPosition(unsigned max);

size_t sumArray(size_t len, int *arr);


void shufflePtrArray(void * data, const size_t len, const size_t sz);





#endif