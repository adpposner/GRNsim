// Reactionrates.c -Calculating reaction rates and times to events (ttes)
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

#include "../include/randomnumbers.h"
#include "../include/globals.h"
#include <string.h>
#include <mkl.h>
#include <float.h>


#define RANDARR_MAX	100

numeric_t_rp logOneMinusU[RND_NUM_BUF_SZ];
const numeric_t_rp * logOneMinusUEnd = logOneMinusU + RND_NUM_BUF_SZ;
numeric_t_rp * currLogOneMinusUPos = logOneMinusU;

VSLStreamStatePtr stream;


#ifdef DOUBLE_PRECISION_RP
 void vLnBuf(numeric_t_rp * srcDest) {
	vdLn(RND_NUM_BUF_SZ,srcDest,srcDest);
}
 void vScalMinus1(numeric_t_rp *srcDest){
	cblas_dscal(RND_NUM_BUF_SZ, -1.0, srcDest, 1);
}
#else
 void vLnBuf(numeric_t_rp * srcDest) {
	vsLn(RND_NUM_BUF_SZ,srcDest,srcDest);
}
 void vScalMinus1(numeric_t_rp *srcDest){
	cblas_sscal(RND_NUM_BUF_SZ, -1.0, srcDest, 1);
}

#endif

void randomUniformDoubleArray(size_t nElems,double min,double max,double *dest){
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, nElems,dest,min,max);
}

void randomUniformSingleArray(size_t nElems,float min,float max, float * dest){
	vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, nElems, dest, min, max);
}


static  void regenLogOneMinusUBuf() {
	
	randomUniformNumberArray(RND_NUM_BUF_SZ, 0.0e0, 0.9999e0, logOneMinusU);
	//The one minus part is really annoying, but I'll stick with it for now
	numeric_t_rp * p = logOneMinusU;
	for(;p!=logOneMinusUEnd;p++)
		*p = 1 - *p;
	vLnBuf(logOneMinusU);
	vScalMinus1(logOneMinusU);
	currLogOneMinusUPos=logOneMinusU;
}

static int randomnumbers_initd = 0;

static void generateRndBufs();

void randomnumbers_init(int seed) {
	assert(seed != 0);
	if(randomnumbers_initd) 
		randomnumbers_free();
	#ifdef FIXSEED
	vslNewStream(&stream, VSL_BRNG_SFMT19937, 42);
	#else
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	#endif
	generateRndBufs();
	randomnumbers_initd = 1;
}


void randomnumbers_free() {
	if(randomnumbers_initd) {
		vslDeleteStream(&stream);
		randomnumbers_initd = 0;
	}
}



void rand_int(size_t lim,int * val) {
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, val, 0, lim);
}

void randomUniqueIntArray(size_t nElems,size_t nMin, size_t nMax,unsigned int *dest) {
	//Using Floyd's algorithm
	unsigned char in_set[RANDARR_MAX] = {0};
	size_t nRange;
	nRange = nMax-nMin;
	if (nRange > RANDARR_MAX) {
		fprintf(stderr,"randomUniqueShortArray: invalid nRange > RANDARRMAX - %zd > %d",nRange,RANDARR_MAX);
		exit(-1);
	} else if (nRange < nElems) {
		fprintf(stderr,"randomUniqueShortArray: range too short for nelems -%zd < %zd",nRange,nElems);
		exit(-1);
	}

	int i,j,currInd;
	currInd = 0;
	for (i=nRange-nElems;i<nRange;i++) {
		viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, 1, &j, 0, i);
		if(in_set[j])
			j = i;
		dest[currInd++]=j+nMin;
		in_set[j]=1;
	}
}

void randomNonUniqueIntArray(size_t nElems,size_t nMin, size_t nMax, unsigned *dest){
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream,nElems,(int *) dest,nMin,nMax);
}


 double getNextExponentialDouble(double rate) {
	if(currLogOneMinusUPos == logOneMinusUEnd){
		regenLogOneMinusUBuf();
	}
	return (*currLogOneMinusUPos++)/rate;
}

 float getNextExponentialSingle(float rate) {
	if(currLogOneMinusUPos == logOneMinusUEnd){
		regenLogOneMinusUBuf();
	}
	return (*currLogOneMinusUPos++)/rate;
}

void randomUniformNumberArray(size_t nElems,effect_t_rp min, effect_t_rp max,effect_t_rp *dest){
#ifdef DOUBLE_PRECISION_RP
randomUniformDoubleArray(nElems, min, max, dest);
#else
randomUniformSingleArray(nElems, min, max, dest);
#endif	
}

rate_t_rp getNextExponential(rate_t_rp rate){
#ifdef DOUBLE_PRECISION_RP
return getNextExponentialDouble(rate);
#else
return getNextExponentialSingle(rate);
#endif	
}

unsigned bindingPositionBuf[RND_NUM_BUF_SZ];
const unsigned * bindingPositionEnd = bindingPositionBuf + RND_NUM_BUF_SZ;
unsigned * currBindingPosition = bindingPositionBuf;



static  void regenBindingPositionBuf(){
	randomNonUniqueIntArray(RND_NUM_BUF_SZ, 0, MAX_NCONNECTIONS, bindingPositionBuf);
	currBindingPosition = bindingPositionBuf;
}


 unsigned getNextBindingPosition(unsigned max) {
	if (currBindingPosition == bindingPositionEnd){
		regenBindingPositionBuf();
	}
	return (*currBindingPosition++) % max;
}


static  void generateRndBufs(){
	regenLogOneMinusUBuf();
	regenBindingPositionBuf();
}

void randomBinomialArray(size_t len, size_t d_low,size_t d_high,double p, unsigned * dest) {
	if (d_high == d_low)
	{
		*dest=d_high;
		return;
	}
	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE, stream, len, (int *) dest, d_high-d_low, p);
	int i;
	for(i=0;i<len;i++)
		dest[i]+=d_low;
}

void randomBimodalBinomialArray(size_t len, size_t d_low, size_t d_high, double p1, double p2, double p_select,
 unsigned * dest) {
	unsigned * dest2;
	unsigned * select_dest;
	if (d_high == d_low)
	{
		*dest=d_high;
		return;
	}
	dest2 = calloc(len, sizeof(*dest2));
	select_dest = calloc(len,sizeof(*select_dest));

	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE, stream, len, (int *) dest, d_high-d_low, p1);
	viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE, stream, len, (int *) dest2, d_high-d_low,p2);
	viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream, len, (int *) select_dest, p_select);
	for(;len;len--){
		dest[len-1] = (select_dest[len-1] ? dest2[len-1] : dest[len-1]) + d_low; 
	}

	free(select_dest);
	free(dest2);
}



size_t sumArray(size_t len, int *arr) {
	unsigned i;
	size_t toRet;
	toRet = 0;
	for (i=0;i<len;i++)
		toRet+=arr[i];
	return toRet;
}



#define SWAP(a,b,size) \
	do {	\
		register size_t __size = (size); \
		register char * __a = (a), *__b = (b);	\
		do {	\
			char __tmp = *__a;	\
			*__a++ = *__b;	\
			*__b++ = __tmp;	\
		}while(--__size > 0); \
	}while(0)	\

void shuffleuCharArray(unsigned char * data, const size_t len){
	int i,j;
	unsigned char tmp;
	for(i=len-1;i>0;i--){
		rand_int(i+1,&j);
		tmp=data[i];
		data[i]=data[j];
		data[j]=tmp;
	}
}

void randomBoolArray(size_t nElems,unsigned int nOnes,unsigned char * dest){
	if (nOnes > nElems){
		fprintf(stderr,"%s:%d - trying to set more miRs active than exist in system, \
			requested %u, setting nactive to nmicro = %zu\n",__FILE__,__LINE__,nOnes,nElems);
		nOnes=nElems;
	}
	memset(dest, 0, nElems);
	memset(dest, 1, nOnes);
	shuffleuCharArray(dest, nElems);
}

void shufflePtrArray(void * data, const size_t len, const size_t sz) {
	char *cdata = (char *) data;
	char *d1,*d2;
	int i,j;
	for(i=len-1;i>0;i--) {
		rand_int(i+1,&j);
		d1 = cdata + (sz * i);
		d2 = cdata + (sz * j);
		SWAP(d1, d2, sz);
	}
}

