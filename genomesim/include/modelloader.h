//modelloader.h - interface for loading model from file
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