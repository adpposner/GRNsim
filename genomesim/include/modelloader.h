#ifndef MODEL_LOADER_H_RP__
#define MODEL_LOADER_H_RP__

#include "globals.h"

extern UnsignedCharArray miRs_states;
struct SimulationComponents;
struct SimOpts;

void model_loader_init(struct SimulationComponents * sim, struct SimOpts *opts);

void loadModel(struct SimulationComponents * sim,struct SimOpts * opts);
void destroyModel(struct SimulationComponents * sim);

#endif