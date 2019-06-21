#ifndef __SIM_H_RP__
#define __SIM_H_RP__
#include "globals.h"
#include "sim/simulator.h"
#include "sim/event.h"
#include "sim/reactionrates.h"
#include "models.h"


enum networkGeneration {USE_FILE,RANDOM_GENERATE};



typedef struct SimulationComponents {
	minHeap *mH;
	GenesList *g;
	ReactionArray *rxns;
} SimulationComponents;


typedef struct SimOpts {

	char jsonFile_[MAX_JSON_FILENAME];
	char dataOutputDir_[MAX_JSON_FILENAME];
    char optSuffix_[MAX_JSON_FILENAME];
	int writeGenes;
	int nactiveMirs;
	int nSimulations;
	unsigned nConformations;
} SimOpts;



#endif
