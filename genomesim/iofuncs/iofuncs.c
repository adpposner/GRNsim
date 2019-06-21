#include "../include/iofuncs/iofuncs.h"
#include "../include/iofuncs/graph_interface.h"
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


#define DEGREESTRINGBUFLEN	5000

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

pmbPtrArray generateDNAPtrArray(GenesList *g){
	pmbPtrArray toRet = pmbPtrArray_alloc(g->nDNA);
	pmbPtrArrayIters it = getpmbPtrArrayIters(&toRet);
	ProducerArrayIters pit = getProducerArrayIters(&g->producers);
	ARRAY_TYPE_FOREACH(pit){
		if(pit.curr->species & DNA)
			*it.curr++ = (pmbPtr){.producer = pit.curr, .species = pit.curr->species};
	}
	assert(it.curr==it.end);
	return toRet;
}

void printDegrees(FILE * degOut, GenesList *g){
    
	pmbPtrArray messes = generateMessengerPtrArray(g);
	pmbPtrArrayIters messIt = getpmbPtrArrayIters(&messes);
	FILE * degreesOut;
	if(!degOut)
		degreesOut = fopen(stProt.degreeFileName, "w");
	else
		degreesOut = degOut;
	//fprintf()
	fprintf(degreesOut,"name\t");
	ARRAY_TYPE_FOREACH(messIt){
		fprintf(degreesOut,"%s\t",messIt.curr->producer->produces.modulator->name);
	}
	fputc('\n', degreesOut);
	fprintf(degreesOut,"mRNA_indegree\t");
	ARRAY_TYPE_FOREACH(messIt){
		fprintf(degreesOut,ltf_rp"\t",getIndegree_mess(messIt.curr->producer));
	}
	fputc('\n', degreesOut);
	char baseProds[DEGREESTRINGBUFLEN];
	char maxProds[DEGREESTRINGBUFLEN];
	char baseDecays[DEGREESTRINGBUFLEN];
	char maxDecays[DEGREESTRINGBUFLEN];
	char * bP = baseProds;
	char * mP = maxProds;
	char * bD = baseDecays;
	char * mD = maxDecays;
	numeric_t_rp baseProdVal, baseDecayVal,maxProdVal,maxDecayVal;
	bP += snprintf(bP,DEGREESTRINGBUFLEN,"mRNA_baseProduction\t");
	mP += snprintf(mP,DEGREESTRINGBUFLEN,"mRNA_maxProduction\t");
	bD += snprintf(bD,DEGREESTRINGBUFLEN,"mRNA_baseDecay\t");
	mD += snprintf(mD,DEGREESTRINGBUFLEN,"mRNA_maxDecay\t");
	ARRAY_TYPE_FOREACH(messIt){
		calcMaxEffect(messIt.curr->producer,&maxProdVal,&baseProdVal,&maxDecayVal,&baseDecayVal);
		bP += sprintf(bP,"%f\t",baseProdVal);
		mP += sprintf(mP,"%f\t",maxProdVal);
		bD += sprintf(bD,"%f\t",baseDecayVal);
		mD += sprintf(mD,"%f\t",maxDecayVal);
	}
	*bP = *mP = *bD = *mD ='\0';
	fprintf(degreesOut,"%s\n",baseProds);
	fprintf(degreesOut,"%s\n",maxProds);
	fprintf(degreesOut,"%s\n",baseDecays);
	fprintf(degreesOut,"%s\n",maxDecays);
	pmbPtrArray genes = generateDNAPtrArray(g);
	pmbPtrArrayIters dnaIt = getpmbPtrArrayIters(&genes);
	memset(baseProds,0,DEGREESTRINGBUFLEN);
	memset(maxProds,0,DEGREESTRINGBUFLEN);
	bP = baseProds;
	mP = maxProds;
	// fprintf(degreesOut,"dna_targ\t");
	// ARRAY_TYPE_FOREACH(dnaIt){
	// 	if(dnaIt.curr->producer->produces.species == MESSENGER){
	// 		fprintf(degreesOut,"%s\t",dnaIt.curr->producer->produces.producer->name);
	// 	}else{
	// 		fprintf(degreesOut,"%s\t",dnaIt.curr->producer->produces.modulator->name);
	// 	}
	// }
	// fputc('\n', degreesOut);
	bP += sprintf(bP,"DNA_baseProduction\t");
	mP += sprintf(mP,"DNA_maxProduction\t");
	ARRAY_TYPE_FOREACH(dnaIt){
		calcMaxEffect(dnaIt.curr->producer,&maxProdVal,&baseProdVal,NULL,NULL);
		bP += sprintf(bP,"%f\t",baseProdVal);
		mP += sprintf(mP,"%f\t",maxProdVal);
	}
	fprintf(degreesOut,"%s\n",baseProds);
	fprintf(degreesOut,"%s\n",maxProds);
	pmbPtrArray_free(&genes);
	pmbPtrArray_free(&messes);
	if(!degOut)
		fclose(degreesOut);
}
struct miniStat fileWriter_init_interval(const char * suffix,species_t species); 
struct miniStat fileWriter_init_unif_time(const char * suffix,species_t species);



void setEventCountInterval(unsigned i){
	stRNA.event_interval=i;
	stProt.event_interval=i;
}

void closeFiles() {
	//fileWriter_close(&stRNA);
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
		fputc('\n', stProt.fh);
 }


 void printReactantQtys_modulator(struct miniStat *st){
	fprintf(st->fh,"%f",clockTime);
		pmbPtrArrayIters it = getpmbPtrArrayIters(&st->data);
		assert(it.start->species & MODULATOR_SPECIES);
			ARRAY_TYPE_FOREACH(it)
				fprintf(st->fh,"\t"ltf_rp,it.curr->modulator->qty);
    fputc('\n',st->fh);
			//fprintf(st->fh,"\t"ltf_rp"\t"ltf_rp"\t"ltf_rp"\n",totalSynths.nMessTranscription,totalSynths.nMicroTranscription,totalSynths.nTranslation);

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




#ifdef UNIFORM_TIME
	 void logData(int k){
		printReactantQtysUnifTime();
	}

	 void setInterval(ulong_type interval){}
	 void printReactantQtysUnifTime(){
	if((clockTime-stProt.lastTime) > stProt.sampling_interval){
	
		printReactantQtys_();
		stProt.lastTime=clockTime;
	}
	}
#elif defined(COMMONEVENTSONLY)
	 void logData(int k){
		printReactantQtysCommonEvents();
	}
	 void setInterval(ulong_type interval){
		stProt.event_interval = interval;
		stRNA.event_interval = interval;
	}
#elif defined(ALLEVENTS)
	 void logData(int k){
		printReactantQtysByInterval(k);
	}
	 void setInterval(ulong_type interval){
		stProt.event_interval = interval;
		stRNA.event_interval = interval;
	}
#else
	#error "Must define UNIFORM_TIME,COMMONEVENTSONLY, or ALLEVENTS"
#endif


void copyConfigfile(const char * srcPath, const char * destPath){
    char srcF[1000];
    char destF[1000];
    sprintf(srcF,"config.xml");
    sprintf(destF,"%s/config.xml",destPath);
    FILE * src, * targ;
    src = fopen(srcF,"r");
    if(!src){
        fprintf(stderr,"cannot open src %s config file\n",srcF);
        exit(12245);
    }
    targ=fopen(destF,"w");
    if(!targ){
        fprintf(stderr,"cannot open targ %s config file\n",destF);
        exit(12245);
    }
    char ch;
    while((ch = fgetc(src))!= EOF)
        fputc(ch,targ);
    fclose(src);
    fclose(targ);
}

void mkdir_if_not_exists(const char * dirName) {
	struct stat st = {0};
	if (stat(dirName, &st) == -1)
		mkdir(dirName, S_IRWXU | S_IRWXG | S_IRWXO | S_IWOTH | S_IROTH);
}

const char * networkParmPath(){
	static char networkPathBuf[1000]= "\0";
	if(!strlen(networkPathBuf)){
		sprintf(networkPathBuf,"output/net.md.%s",networkParamHash(0));
		mkdir_if_not_exists(networkPathBuf);}
	return networkPathBuf;
}

const char * specificNetworkPath(){
	static char specificPathBuf[1000]="\0";
	if(!strlen(specificPathBuf)){
		sprintf(specificPathBuf, "%s/ngs.%d",networkParmPath(),networkGenSeed);
		mkdir_if_not_exists(specificPathBuf);}
	return specificPathBuf;
}

const char * networkPathWNActiveMiRs(unsigned nActiveMirs){
	static char networkPathWNActiveMiRsBuf[1000]="\0";
	if (!strlen(networkPathWNActiveMiRsBuf)){
		sprintf(networkPathWNActiveMiRsBuf, "%s/nMiRs.%zd",specificNetworkPath(),nActiveMirs);
		mkdir_if_not_exists(networkPathWNActiveMiRsBuf);}
	return networkPathWNActiveMiRsBuf;
}
static char jsonfilepth[1000];

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

char * writeGenesList(GenesList * g,SkipDisabled skipDisabled,char * root,unsigned nActive,char * miRsStates) {
	static char jsonfile[1000];
	int written;
	if(skipDisabled == SKIP_DISABLED){
		written = snprintf(jsonfilepth, 1000, "%s/miRs.%d",root,nActive);
		mkdir_if_not_exists(jsonfilepth);
		written += snprintf(jsonfilepth+written,1000,"/states.%s",miRsStates);
		mkdir_if_not_exists(jsonfilepth);

		snprintf(jsonfile,1000,"%s/genes.json",jsonfilepth);

           
	} else {
		snprintf(jsonfile,1000,"%s/genes.json",specificNetworkPath());
		
	}
    copyConfigfile(".",specificNetworkPath());
	writeGenesList_JSON(g, jsonfile,skipDisabled);
	return jsonfile;
}

GenesList * readGenesList(const char * filename){
    static char readGenesListBuf[1000];
    char * fname, * dname;
    
    //split_file_from_path(filename,&dname,&fname);
    if(!strstr(filename,".json")){
        snprintf(readGenesListBuf,1000,"%s/genes.json",filename);
    }else{
        snprintf(readGenesListBuf,1000,"%s",filename);
    }
	return readGenesList_JSON(readGenesListBuf);
}

char * outputFileString(struct miniStat * st) {

	sprintf(st->dirName,"%s",jsonfilepth);
	
	sprintf(st->fileName,"%s/sim.%d.%s.txt",st->dirName,simulationSeed,st->suffix);
	sprintf(st->degreeFileName,"%s/sim.%d.%s.degrees.txt",
		st->dirName,simulationSeed,st->suffix);
	return st->fileName;
}






// create JSON output filename string
char * JSONOutputFileString(struct miniStat *st,UnsignedCharArray miRs_active){
	static char strOut[1000];
	sprintf(strOut,"%s/sim.%d.mirs.%d.json",st->dirName,simulationSeed,UnsignedCharArray_total(miRs_active));
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
struct miniStat fileWriter_init_unif_time(const char * suffix,species_t species) {
	struct miniStat toRet = fileWriter_init_by_unif_time(suffix, species);
	outputFileString(&toRet);
	mkdir_if_not_exists(toRet.dirName);
	fileWriter_open_(&toRet);
	return toRet;
}




void openFilesForWriting(GenesList *g,char * optSuffix) {
	stRNA = fileWriter_init_unif_time("rna",MESSENGER);
	stProt = fileWriter_init_unif_time(optSuffix ? optSuffix : "prot",PROTEIN);
	stRNA.data = generateMessengerPtrArray(g);
	stProt.data = generateProteinPtrArray(g);
	printDegrees(NULL,g);
}


#endif
