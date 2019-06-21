#ifndef MODEL_GENERATOR_H__
#define MODEL_GENERATOR_H__
#include "models.h"
#include "sim.h"

void createNetwork(SimulationComponents *sim);
void destroyNetwork(SimulationComponents *toFree);
#endif