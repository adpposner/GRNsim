//iofuncs.c - source for file I/O for both netgen and sim
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
#include "../include/iofuncs/iofuncs.h"
#include "../include/iofuncs/json_interface.h"
#include "../include/globals.h"
#include "../include/models/geneslist_p.h"
#include "../include/models/basicgeneticelement_p.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../include/globals.h"

#define XSTR(name) #name
#define DEFINE_SIMPLE_ARRAY(name) name##Array name ##Array_alloc(size_t len) \
{ name##Array toRet; toRet.data = malloc(len*sizeof(*toRet.data)); toRet.len = len; return toRet;}	\
void name ##Array_free(name##Array *toFree) {if(toFree->data) \
{free(toFree->data);toFree->data = NULL;}toFree->len=0;} \

#define DEFINE_WRITE_FUNCS(name,fmtChar) int name##Array_write_binary(struct miniStat *st,name##Array arr){ \
fwrite(&arr.len,sizeof(arr.len),1,st->fh); \
return fwrite(arr.data,sizeof(*arr.data),arr.len,st->fh);} \
void name##Array_print(struct miniStat *st,name##Array arr) { \
char buf[5000]; char *pch = buf; \
int i; for(i=0;i<arr.len; i++) { \
	pch += sprintf(pch,"%"fmtChar"\t",arr.data[i]); \
}*--pch = '\n';fputs(buf, st->fh); }\


struct miniStat stRNA, stProt,stOutput;
int idxLast=0;




pmbPtrArray generateMessengerPtrArray(GenesList *g){
	pmbPtrArray toRet = pmbPtrArray_alloc(g->nMess);
	pmbPtrArrayIters it = getpmbPtrArrayIters(&toRet);
	ProducerArrayIters pit = getProducerArrayIters(&g->producers);
	ARRAY_TYPE_FOREACH(pit){
		if(pit.curr->species & MESSENGER)
			*it.curr++ = (pmbPtr){.producer = pit.curr, .species = pit.curr->species};
	}
	assert(it.curr==it.end);
	return toRet;
}

pmbPtrArray generateProteinPtrArray(GenesList *g){
	pmbPtrArray toRet = pmbPtrArray_alloc(g->nProt);
	pmbPtrArrayIters it = getpmbPtrArrayIters(&toRet);
	ModulatorArrayIters mit = getModulatorArrayIters(&g->modulators);
	ARRAY_TYPE_FOREACH(mit){
		if(mit.curr->species & PROTEIN)
			*it.curr++ = (pmbPtr){.modulator = mit.curr, .species = mit.curr->species};
	}
	assert(it.curr==it.end);
	return toRet;
}

void printDegrees(GenesList *g){
	pmbPtrArray messes = generateMessengerPtrArray(g);
	pmbPtrArrayIters messIt = getpmbPtrArrayIters(&messes);
	FILE * degreesOut;
	degreesOut = fopen(stProt.degreeFileName, "w");
	ARRAY_TYPE_FOREACH(messIt){
		fprintf(degreesOut,"%s\t",messIt.curr->producer->produces.modulator->name);
	}
	fputc('\n', degreesOut);
	ARRAY_TYPE_FOREACH(messIt){
		fprintf(degreesOut,ltf_rp"\t",getIndegree_mess(messIt.curr->producer));
	}
	fputc('\n',degreesOut);
	pmbPtrArray_free(&messes);
	fclose(degreesOut);
}

struct miniStat fileWriter_init_interval(const char * suffix,species_t species); 
struct miniStat fileWriter_init_unif_time(const char * suffix,const char * outDir, species_t species);



void setEventCountInterval(unsigned i){
	stRNA.event_interval=i;
	stProt.event_interval=i;
}

void closeFiles() {
	fileWriter_close(&stProt);
	fileWriter_close(&stRNA);
}



 void printHeaders(GenesList *g){
		GenesListIterator gI = getGenesListIters(g);		  	 
		fprintf(stProt.fh,"%s","clockTime");

		ARRAY_TYPE_FOREACH(gI.mIt){
			if(gI.mIt.curr->species & (PROTEIN))
				fprintf(stProt.fh,"\t%s",gI.mIt.curr->name);
		}
        fprintf(stProt.fh,"\ttotal_mess_txc\ttotal_micro_txc\ttotal_txl\n");
 }


 void printReactantQtys_modulator(struct miniStat *st){
	fprintf(st->fh,"%f",clockTime);
		pmbPtrArrayIters it = getpmbPtrArrayIters(&st->data);
		assert(it.start->species & MODULATOR_SPECIES);
			ARRAY_TYPE_FOREACH(it)
				fprintf(st->fh,"\t"ltf_rp,it.curr->modulator->qty);
			fprintf(st->fh,"\t"ltf_rp"\t"ltf_rp"\t"ltf_rp"\n",totalSynths.nMessTranscription,
            totalSynths.nMicroTranscription,totalSynths.nTranslation);

}

 void printReactantQtys_producer(struct miniStat *st){
	fprintf(st->fh,"%f",clockTime);
		pmbPtrArrayIters it = getpmbPtrArrayIters(&st->data);
		assert(it.start->species & PRODUCER_SPECIES);
			ARRAY_TYPE_FOREACH(it)
				fprintf(st->fh,"\t"ltf_rp,it.curr->producer->qty);
			fputc('\n', st->fh);
}




 void printReactantQtys_(){
#ifdef PRINT_RNA
     printReactantQtys_producer(&stRNA);
#endif
	printReactantQtys_modulator(&stProt);
}

 void logStart(){
	printReactantQtys_();
}

 void printReactantQtysByInterval(int idx) {
	 if ((idx % stProt.event_interval) == 0)
			printReactantQtys_();
		
}

 void printReactantQtysCommonEvents() {
 	if (commonBirthDeaths == commonBirthDeathsLast){
 		return;
 	}else if (((commonBirthDeaths-commonBirthDeathsLast) % stProt.event_interval) == 0){
		printReactantQtys_();
		commonBirthDeathsLast = commonBirthDeaths;
	}
}




	 void logData(int k){
		printReactantQtysUnifTime();
	}


	 void printReactantQtysUnifTime(){
	if((clockTime-stProt.lastTime) > stProt.sampling_interval){
	
		printReactantQtys_();
		stProt.lastTime=clockTime;
	}
	}


void mkdir_if_not_exists(const char * dirName) {
	struct stat st = {0};
	if (stat(dirName, &st) == -1)
		mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO | S_IWOTH | S_IROTH);
}


static void split_file_from_path(char * fnamestring, char ** dirname, char ** fname){
    static char localdirname[MAX_JSON_FILENAME];
    static char localfilename[MAX_JSON_FILENAME];
    *localdirname = *localfilename = '\0';
    char * sep = fnamestring;
    char * nextItem;
    while ((nextItem = strpbrk(sep+1,"\\/"))) sep = nextItem;
    if (fnamestring != sep) sep++;
    strncat(localdirname,fnamestring,sep-fnamestring);
    strncat(localfilename,sep,MAX_JSON_FILENAME);
   // printf("%s:%d\t filename = %s\t dirname = %s\n",__FILE__,__LINE__,localfilename,localdirname);
    *dirname = localdirname;
    *fname = localfilename;
}

//NOT USED FOR MAIN PAPER - ONLY FOR DEMO

char * writeGenesList(GenesList * g,SkipDisabled skipDisabled, const char * dest) {
	static char jsonfile[1000];
    snprintf(jsonfile,1000,"%s",dest);

	writeGenesList_JSON(g, jsonfile,skipDisabled);
	return jsonfile;
}

char  * writeTempGenesList(GenesList * g, const char * infile,const char * destdir){
    static char tmpoutfilename[1000];
    char * inputDir, *inputFilename;
    split_file_from_path(infile,&inputDir,&inputFilename);
    snprintf(tmpoutfilename,MAX_JSON_FILENAME,"%s/%s",destdir,inputFilename);
    mkdir_if_not_exists(destdir);
    writeGenesList_JSON(g,tmpoutfilename,SKIP_DISABLED);
    return tmpoutfilename;
}


GenesList * readGenesList(const char * filename){
	return readGenesList_JSON(filename);
}

char * outputFileString(struct miniStat * st,const char * outdir) {
    int status;
	status = snprintf(st->dirName,MAX_JSON_FILENAME,"%s",outdir);
	if(status < 0){
        fprintf(stderr,"Failure in func snprintf at %s:%d\n",__FILE__,__LINE__);
    }
	status = snprintf(st->fileName,MAX_JSON_FILENAME,"%s/sim.%d.%s.txt",st->dirName,simulationSeed,st->suffix);
	if(status < 0){
        fprintf(stderr,"Failure in func snprintf at %s:%d\n",__FILE__,__LINE__);
    }
    status = snprintf(st->degreeFileName,MAX_JSON_FILENAME,"%s/sim.%d.%s.degrees.txt",
		st->dirName,simulationSeed,st->suffix);
	if(status < 0){
        fprintf(stderr,"Failure in func snprintf at %s:%d\n",__FILE__,__LINE__);
    }
    return st->fileName;
}

// create JSON output filename string
char * JSONOutputFileString(struct miniStat *st,UnsignedCharArray miRs_active){
	static char strOut[1000];
	snprintf(strOut,1000,"%s/sim.%d.mirs.%d.json",st->dirName,simulationSeed,UnsignedCharArray_total(miRs_active));
	return strOut;
}

#include "../include/modelloader.h"
void saveStateJSON(GenesList *g,SkipDisabled skipDisabled){
	writeGenesList_JSON(g, JSONOutputFileString(&stProt,miRs_states),skipDisabled);
}

void fileWriter_open_(struct miniStat * st) {
	st->fh = fopen(st->fileName,"wb");
}





void fileWriter_close(struct miniStat *st) {
	fclose(st->fh);
	st->fh = NULL;
	pmbPtrArray_free(&st->data);
	st->data = emptypmbPtrArray();
}
#if defined(COMMONEVENTSONLY) || defined(ALLEVENTS)
struct miniStat fileWriter_init_by_event_count(const char * suffix,species_t species){
	struct miniStat toRet = {.fh = NULL, .event_interval = INITLOGGINGINTERVAL,.species=species,.lastTime=0.0};
	strncpy(toRet.suffix, suffix, strlen(suffix));
	return toRet;
}

struct miniStat fileWriter_init_interval(const char * suffix,species_t species) {
	struct miniStat toRet = fileWriter_init_by_event_count(suffix,species);
	outputFileString(&toRet);
	mkdir_if_not_exists(toRet.dirName);
	fileWriter_open_(&toRet);
	return toRet;
}



void openFilesForWriting(GenesList *g,char * optSuffix) {
	stRNA = fileWriter_init_interval("rna",MESSENGER);
	stProt = fileWriter_init_interval(optSuffix ? optSuffix : "prot",PROTEIN);
	stRNA.data = generateMessengerPtrArray(g);
	stProt.data = generateProteinPtrArray(g);
	totalSynths.nMessTranscription = totalSynths.nMicroTranscription = totalSynths.nTranslation = 0;
	printDegrees(g);
}



#else

struct miniStat fileWriter_init_by_unif_time(const char * suffix,species_t species){
	struct miniStat toRet = {.fh = NULL,.sampling_interval=UNIFORMSAMPLINGINTERVAL ,.species=species,.lastTime=0.0};
	strncpy(toRet.suffix, suffix, strlen(suffix));
	return toRet;
}
struct miniStat fileWriter_init_unif_time(const char * suffix,const char * outDir, species_t species) {
	struct miniStat toRet = fileWriter_init_by_unif_time(suffix, species);
	outputFileString(&toRet,outDir);
	mkdir_if_not_exists(toRet.dirName);
	fileWriter_open_(&toRet);
	return toRet;
}


void openFilesForWriting(GenesList *g,const char * destDir,char * optSuffix) {
    mkdir_if_not_exists(destDir);
	stRNA = fileWriter_init_unif_time("rna",destDir,MESSENGER);
	stProt = fileWriter_init_unif_time(optSuffix ? optSuffix : "prot",destDir,PROTEIN);
	stRNA.data = generateMessengerPtrArray(g);
	stProt.data = generateProteinPtrArray(g);
	printDegrees(g);
}


#endif
