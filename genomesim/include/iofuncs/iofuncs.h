#ifndef IOFUNCS_H_RP__
#define IOFUNCS_H_RP__
#include "../models/geneslist.h"
#include "../parameters.h"


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
	char dirName[256];
	char fileName[256];
	char suffix[56];
	char degreeFileName[256];
	union {
		ulong_type event_interval;
		time_t_rp sampling_interval;
	};
	pmbPtrArray data;
	species_t species;
	time_t_rp lastTime;
};

extern struct miniStat stRNA, stProt;


void openFilesForWriting(GenesList *g,char * optSuffix);
void setInterval(ulong_type i);
void closeFiles();

// void printReactantQtys( struct miniStat *st,int idx);
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
void printDegrees(FILE * degOut, GenesList *g);
char * writeGenesList(GenesList * g,SkipDisabled skipDisabled,char * root,unsigned nActive,char * miRsStates);
//void writeGraph(connectedGenesList * cg);

GenesList * readGenesList(const char * filename);


#endif