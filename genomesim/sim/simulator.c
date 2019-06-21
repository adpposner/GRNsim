
// simulator.c - source for interface to mainsim
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
#include "../include/sim.h"
#include "../include/iofuncs.h"
#include "../include/randomnumbers.h"
#include <string.h>

#ifdef DEBUG 
#define INITNCYCLES  10000
#ifndef FIXED_N_CYCLES
#define FIXED_N_CYCLES	100000
#define RUN_UNTIL	100.0
#endif
#elif defined(FEW_CYCLES)
#define INITNCYCLES  10000
#ifndef FIXED_N_CYCLES
#define FIXED_N_CYCLES	100000
#define RUN_UNTIL	100.0
#endif
#else
#define INITNCYCLES  10000
#define FIXED_N_CYCLES	75000000
#ifndef RUN_UNTIL
#define RUN_UNTIL	10000.0
#endif
#endif



void setupSimulation(SimulationComponents *sim,SimOpts * opts) {
	time_t t;
	srand((unsigned) time(&t));
	#ifdef DEBUG
	simulationSeed = 252;
	#else
	simulationSeed = rand() % 10000000;
	#endif
	randomnumbers_init(simulationSeed);
	resetTime();
	resetInitialQuantities(sim->g);
	initAllRatesTimes(sim->rxns, sim->mH);
	testMinHeap(*sim->mH);
    printf("dataoutputdir=%s\n",opts->dataOutputDir_);
	openFilesForWriting( sim->g, strlen(opts->dataOutputDir_) ? opts->dataOutputDir_ : ".",
    strlen(opts->optSuffix_) ? opts->optSuffix_ : NULL);
	printHeaders(sim->g);
	logStart();
}


void shutDownSimulation(SimulationComponents *sim) {
	closeFiles();
	randomnumbers_free();
}

// static  void runSim(SimulationComponents* sim) {
// 	int k;
// 	ulong_type initialDuration = INITNCYCLES;
// 	ulong_type loggingInterval = DEFAULTLOGGINGINTERVAL;

// 	setInterval(INITLOGGINGINTERVAL);
// 	for (k = 0; k < initialDuration; k++)
// 	{
// 		executeReaction(sim->mH); 
// 		logData(k);}


// 	while (commonBirthDeaths < FIXED_N_CYCLES)
// 	{executeReaction(sim->mH); logData(k);}
// }

 void runSimUntil(SimulationComponents *sim, double until) {
	unsigned k = 0;
	while (clockTime < until)
	{executeReaction(sim->mH); logData(k); k++;}
}

 void doSimulation(SimulationComponents * sim,SimOpts * opts) {
     int nSimulations = opts->nSimulations;
	for (;nSimulations > 0;nSimulations--)
	{	
		setupSimulation(sim,opts);
		runSimUntil(sim, RUN_UNTIL);
		shutDownSimulation(sim);
		printTotals();
	}

}



