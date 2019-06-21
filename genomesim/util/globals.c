// globals.c - source for defs of global/debug functions and variables
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
#include "../include/globals.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


#ifdef DEBUG
char * species_names(species_t sp) {
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
#endif

void UnsignedInt_free(unsigned int * toFree){}
void UnsignedChar_free(unsigned char *toFree){}

DEFINEBASICARRAYTYPEPRIMITIVE(unsigned int,UnsignedInt)
DEFINEBASICARRAYTYPEPRIMITIVE(unsigned char,UnsignedChar)
ulong_type mallocDBGCt=0;
ulong_type freeDBGCt=0;


int simulationSeed = -1;
int miRNAs_active = -1;
int networkGenSeed = -1;

time_t_rp clockTime=0.0;

void resetTime(){clockTime=0.0;commonBirthDeaths=0;
	totalSynths.nMessTranscription = totalSynths.nMicroTranscription =totalSynths.nTranslation = 0;}

void * mallocDB(size_t sz,const char * caller)	{
mallocDBGCt++;return malloc(sz);}
void freeDB(void * toFree,const char * caller)		{
free(toFree);freeDBGCt++;}


int Unsignedcmpfunc(const void *a, const void *b) {
	unsigned *l,*r;
	l=(unsigned*)a;r=(unsigned*)b;
	if (*l > *r) return 1;
	else if (*r > *l) return -1;
	else return 0;
}

void UnsignedIntArray_sort(UnsignedIntArray arr) {
	qsort(arr.data, arr.length, sizeof(*arr.data), Unsignedcmpfunc);
}

UnsignedIntArray UnsignedIntArray_fromArray(unsigned * start, unsigned *end){
	size_t sz = end-start;
	UnsignedIntArray toRet;
	printf("%s\n",__func__);
	toRet = UnsignedIntArray_alloc(sz);
	memcpy(toRet.data, start, sz * sizeof(*toRet.data));
	return toRet;
}

UnsignedIntArray UnsignedIntArray_cat(UnsignedIntArray arr1, UnsignedIntArray arr2){
	UnsignedIntArray toRet;
	ulong_type totalLength;
	totalLength = arr1.length + arr2.length;
	printf("%s\n",__func__);
	toRet = UnsignedIntArray_alloc(totalLength);
	memcpy(toRet.data, arr1.data, arr1.length * sizeof(*arr1.data));
	memcpy(toRet.data + arr1.length, arr2.data, arr2.length * sizeof(*arr2.data));
	return toRet;
}

unsigned UnsignedIntArray_total(UnsignedIntArray arr) {
	unsigned x;
	x=0;
	unsigned i;
	for(i=0;i<arr.length;i++)
		x += arr.data[i];
	return x;
}

void UnsignedCharArray_clear(UnsignedCharArray arr){
	memset(arr.data, 0, arr.length*sizeof(*arr.data));
}

unsigned UnsignedCharArray_total(UnsignedCharArray arr){
	unsigned toRet;
	toRet=0;
	unsigned i;
	for(i=0;i<arr.length;i++)
		toRet += arr.data[i] ? 1 : 0;
	return toRet;
}

char * readFileToString(const char * filename, size_t * strSize) {
	char *strIn;
	FILE * inFile;
	inFile = fopen(filename,"rb");
	fseek(inFile,0L,SEEK_END);
	*strSize = ftell(inFile);
	fseek(inFile, 0L, SEEK_SET);
	strIn = calloc(sizeof(*strIn),*strSize+1);
	fread(strIn,sizeof(*strIn),*strSize,inFile);
	fclose(inFile);
	return strIn;
}

