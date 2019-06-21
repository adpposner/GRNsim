
#include "../include/sim.h"
#include "../include/iofuncs.h"
#include "../include/randomnumbers.h"


#ifdef DEBUG 
#define INITNCYCLES  10000
#ifndef FIXED_N_CYCLES
#define FIXED_N_CYCLES	100000
#define RUN_UNTIL	100.0
#endif
//#define NCYCLES 2000
#elif defined(FEW_CYCLES)
#define INITNCYCLES  10000
#ifndef FIXED_N_CYCLES
#define FIXED_N_CYCLES	100000
#define RUN_UNTIL	100.0
#endif
#else
#define INITNCYCLES  10000
#define FIXED_N_CYCLES	275000000
#ifndef RUN_UNTIL
#define RUN_UNTIL	2500.0
#endif
#endif



void setupSimulation(SimulationComponents *sim) {
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
	char * suffix = NULL;
	openFilesForWriting( sim->g, suffix);
	printHeaders(sim->g);
	logStart();
}


void shutDownSimulation(SimulationComponents *sim) {
	closeFiles();
	randomnumbers_free();
}

static  void runSim(SimulationComponents* sim) {
	int k;
	ulong_type initialDuration = INITNCYCLES;
	ulong_type loggingInterval = DEFAULTLOGGINGINTERVAL;

	setInterval(INITLOGGINGINTERVAL);
	for (k = 0; k < initialDuration; k++)
	{
		executeReaction(sim->mH); 
		logData(k);}
		setInterval(loggingInterval);

	while (commonBirthDeaths < FIXED_N_CYCLES)
	{executeReaction(sim->mH); logData(k);}
}

 void runSimUntil(SimulationComponents *sim, double until) {
	unsigned k = 0;
	while (clockTime < until)
	{executeReaction(sim->mH); //logData(k); k++;
	executeReaction(sim->mH); //logData(k); k++;
    executeReaction(sim->mH);
    executeReaction(sim->mH);
    executeReaction(sim->mH);
    executeReaction(sim->mH);
    logData(k);k+=8;}
}

 void doSimulation(SimulationComponents * sim,int nSimulations) {
	for (;nSimulations > 0;nSimulations--)
	{	
		setupSimulation(sim);
		runSimUntil(sim, RUN_UNTIL);
		shutDownSimulation(sim);
		printTotals();
	}

}



