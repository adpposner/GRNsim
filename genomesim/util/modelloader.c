// modelloader.c - source for interface between model generation and model simulation
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
#include "../include/modelloader.h"
#include "../include/models.h"
#include "../include/randomnumbers.h"
#include "../include/iofuncs/iofuncs.h"
#include "../include/sim.h"
#include "../include/models/basicgeneticelement_p.h"


UnsignedCharArray miRs_states;

void printmiRsStats(){
	UnsignedCharArrayIters it = getUnsignedCharArrayIters(&miRs_states);
	ARRAY_TYPE_FOREACH(it) printf("miR is %s\n",(*it.curr) ? "ENABLED" : "DISABLED");
}

void init_miRs_states(unsigned nMicro){
	miRs_states = UnsignedCharArray_alloc(nMicro);
	UnsignedCharArray_clear(miRs_states);
}
char * miRsStatesToString(){
	static char buf[500] = {0};
	UnsignedCharArrayIters it = getUnsignedCharArrayIters(&miRs_states);
	char * pch = buf;
	ARRAY_TYPE_FOREACH(it) *pch++ = *it.curr + 'A';
	*pch='\0';
	return buf;
}


static void toggleMir(Modulator * micro,unsigned char activate){
	assert(micro->species == MICRO);
	if(activate){
		micro->isEnabled = 1;
	} else {
		micro->isEnabled = 0;
	}
}

void toggleFamilyForMirs(GenesList * g){

}

void setMiRNAsFromStates(GenesList *g){
	pmbPtrArray microArr = getPtrArrayForType(g, MICRO);
	if(microArr.length != miRs_states.length) 
		{fprintf(stderr,"Mismatch in mir array lengths %zu != %zu",microArr.length,miRs_states.length);exit(-782);}
	pmbPtrArrayIters mIt =getpmbPtrArrayIters(&microArr);
	UnsignedCharArrayIters stateIt = getUnsignedCharArrayIters(&miRs_states);
	for(mIt.curr=mIt.start,stateIt.curr=stateIt.start;mIt.curr!=mIt.end;mIt.curr++,stateIt.curr++){
		toggleMir(mIt.curr->modulator, *stateIt.curr);
	}
	assert(stateIt.curr==stateIt.end);
	pmbPtrArray_free(&microArr);
	toggleFamilyForDisabled(g);
}


void setRandomMiRNAsActive(GenesList *g,ulong_type nActive){
	randomBoolArray(miRs_states.length, nActive, miRs_states.data);
	setMiRNAsFromStates(g);
}


void disableMiRNAs(GenesList *g){
	UnsignedCharArray_clear(miRs_states);
	setMiRNAsFromStates(g);
}

void enableAllMiRNAs(GenesList *g){
	UnsignedCharArrayIters it = getUnsignedCharArrayIters(&miRs_states);
	ARRAY_TYPE_FOREACH(it) *it.curr=1;
	setMiRNAsFromStates(g);
}

void model_loader_init(SimulationComponents * sim,SimOpts * opts){
	ulong_type totalNMicro = nMicro_total(sim->g);
	totalSynths.nTranslation = totalSynths.nMessTranscription = totalSynths.nMicroTranscription = 0;
	init_miRs_states(totalNMicro);
	if (opts->nactiveMirs == 0)
		disableMiRNAs(sim->g);
	else if (opts->nactiveMirs == totalNMicro)
		enableAllMiRNAs(sim->g);
	else
		setRandomMiRNAsActive(sim->g, opts->nactiveMirs);
    char  * tmpfilepath;
    tmpfilepath = writeTempGenesList(sim->g,opts->jsonFile_,opts->dataOutputDir_);
	//char * filePath = writeGenesList(sim->g, SKIP_DISABLED, opts->jsonFile_);
   // fprintf(stdout,"%s:%d\t tmp file path = %s\n",__FILE__,__LINE__,tmpfilepath);
	genesList_free(sim->g);
	sim->g = readGenesList(tmpfilepath);	
}



#include "../include/iofuncs.h"
void loadModel(SimulationComponents * sim,SimOpts * opts){
	time_t t;
	srand((unsigned) time(&t));
	unsigned activeMirSeed = rand() % 10000000;
	randomnumbers_init(activeMirSeed);
	char filePath[1000];
	int status;
    status = snprintf(filePath,1000,"%s",opts->jsonFile_);
    if(status < 0){
        fprintf(stderr,"Failure in func snprintf at %s:%d\n",__FILE__,__LINE__);
    }
	sim->g = readGenesList(filePath);
	model_loader_init(sim, opts);
	generateReactionsAndDependencies(sim->g, &sim->rxns);
	sim->mH = minHeap_init(sim->rxns);
	randomnumbers_free();
}



void destroyModel(SimulationComponents * sim){
	minHeap_free(sim->mH);
	free(sim->mH);
	sim->mH = NULL;
	ReactionArray_free(sim->rxns);
	free(sim->rxns);
	sim->rxns = NULL;
	UnsignedCharArray_free(&miRs_states);
	genesList_free(sim->g);
	sim->g = NULL;
}
