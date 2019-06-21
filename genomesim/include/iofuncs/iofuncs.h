//iofuncs.h - header for file I/O for both netgen and sim
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
#ifndef IOFUNCS_H_RP__
#define IOFUNCS_H_RP__
#include "../models/geneslist.h"
#include "../parameters.h"
#include "../globals.h"


#ifndef MAX_JSON_FILENAME
#define MAX_JSON_FILENAME 500
#endif


#define CPU_TIMEIT_START	clock_t begin; begin=clock();
#define CPU_TIMEIT_END		clock_t end=clock(); \
printf("time spent: %f\n",(double)(end-begin)/CLOCKS_PER_SEC);
#define WALL_TIMEIT_START	time_t beg = time(NULL);
#define WALL_TIMEIT_END		printf("wall time spend: %f\n",(double)time(NULL)-beg);

#define ARRAY_TYPE_NAME(name)	name ## Array
#define DECLARE_SIMPLE_ARRAY(name)	typedef struct name ##Array { \
	name * data;	\
	unsigned len;	\
} name ## Array;	\
\
name##Array name##Array_alloc(size_t len); \
void name##Array_free(name##Array *toFree);\

#define DECLARE_WRITE_FUNCS(name) int name##Array_write_binary(struct miniStat *st, name##Array arr); \
void name##Array_print(struct miniStat *st,name##Array arr); \

struct miniStat {
	FILE * fh;
	char dirName[MAX_JSON_FILENAME];
	char fileName[MAX_JSON_FILENAME];
	char suffix[MAX_JSON_FILENAME];
	char degreeFileName[MAX_JSON_FILENAME];
	union {
		ulong_type event_interval;
		time_t_rp sampling_interval;
	};
	pmbPtrArray data;
	species_t species;
	time_t_rp lastTime;
};

extern struct miniStat stRNA, stProt;


void openFilesForWriting(GenesList *g,const char * destDir,char * optSuffix);

void closeFiles();


 void printReactantQtysProtEvents( struct miniStat *st);
 void writeReactantQtys(struct miniStat *st,int idx, GenesList g);
 void printReactantQtysByInterval(int idx);
 void printReactantQtysCommonEvents(); 
 void printReactantQtysUnifTime();

DECLARE_SIMPLE_ARRAY(double)
DECLARE_WRITE_FUNCS(double)

DECLARE_SIMPLE_ARRAY(int)
DECLARE_WRITE_FUNCS(int)

struct miniStat fileWriter_init(const char * suffix, const ulong_type loggingInterval,species_t species);
void fileWriter_close(struct miniStat *st);

void printHeaders(GenesList *g);



char * writeGenesList(GenesList * g,SkipDisabled skipDisabled,const char * dest);
char  * writeTempGenesList(GenesList * g, const char * infile,const char * destdir);


GenesList * readGenesList(const char * filename);


#endif