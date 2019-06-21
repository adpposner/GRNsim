//simulator.h - Interface for running sim
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
#ifndef SIMULATOR_H_RP__
#define SIMULATOR_H_RP__

struct SimulationComponents;
struct SimOpts;
void setupSimulation(struct SimulationComponents *sim,struct SimOpts * opts);
void shutDownSimulation(struct SimulationComponents *sim);
 void doSimulation(struct SimulationComponents *sim, struct SimOpts * opts);

#endif