//modelgenerator.h - interface for file I/O for both netgen and sim
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
#ifndef MODEL_GENERATOR_H__
#define MODEL_GENERATOR_H__
#include "models.h"
#include "sim.h"

void createNetwork(SimulationComponents *sim,const char * dest);
void destroyNetwork(SimulationComponents *toFree);
#endif