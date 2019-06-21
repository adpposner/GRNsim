//sim.h - Primary structure for simulation options
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
#ifndef __SIM_H_RP__
#define __SIM_H_RP__


#include "sim/simulator.h"
#include "sim/event.h"
#include "sim/reactionrates.h"
#include "models.h"
#include "globals.h"



enum networkGeneration {USE_FILE,RANDOM_GENERATE};

#ifndef MAX_JSON_FILENAME
#define MAX_JSON_FILENAME 500
#endif

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