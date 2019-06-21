//geneslist.h -Structures/functions for geneslist
//The geneslist is the primary data structure for the simulation. 
//It contains network info and parameters, but the only "dynamic" value is quantity
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
#ifndef GENESLIST_RP_H__
#define GENESLIST_RP_H__

#include "../globals.h"
#include "basicgeneticelement.h"




typedef struct GenesList GenesList;


//resets ICs for new sim
void resetInitialQuantities(GenesList *g);

pmbPtrArray getPtrArrayForType(GenesList *g, species_t sptype);

GenesList * genesList_alloc(const ulong_type nMess, const ulong_type nMicro, 
	const ulong_type nMessMir, const ulong_type nTFDNA,unsigned char allocInternals);

GenesList * genesList_base_generate(ulong_type nMess, ulong_type nMicro);
void genesList_bound_quantities_init(GenesList *g, const ulong_type nMessMir, const ulong_type nTFDNA);
void assignBoundElements(GenesList *g);
void initDefaultGenesListQuantities(GenesList *g);

void genesList_free(GenesList * tofree);


ulong_type nMicro_total(GenesList * g);

ulong_type nEnabledForType(GenesList * g,species_t species);

void toggleFamilyForDisabled(GenesList *g);
unsigned isEnabled(const pmbPtr pmb);
unsigned char isProducerEnabled(Producer * producer);
unsigned char isModulatorEnabled(Modulator * modulator);
unsigned char isBoundEnabled(BoundElement * bd);
#endif