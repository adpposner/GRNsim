
#ifndef SIMULATOR_H_RP__
#define SIMULATOR_H_RP__

struct SimulationComponents;
struct SimOpts;
void setupSimulation(struct SimulationComponents *sim);
void shutDownSimulation(struct SimulationComponents *sim);
 void doSimulation(struct SimulationComponents *sim,int nSimulations);

#endif