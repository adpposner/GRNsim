//geneslist_p.h - "Private" structures/functions for geneslist
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
#ifndef GENESLIST_P_H__
#define GENESLIST_P_H__

#include "basicgeneticelement.h"

typedef struct initialQuantities{
	UnsignedIntArray prodICs;
	UnsignedIntArray modICs;
	UnsignedIntArray boundICs;
} initialQuantities;

typedef struct GenesList {
	ProducerArray producers;
	ModulatorArray modulators;
	BoundElementArray bounds;
	ulong_type nDNA, nMess, nMicro,nProt,nMessMir,nTFDNA;
	initialQuantities * ICs;
} GenesList;


typedef struct GenesListIterator {
	ProducerArrayIters pIt;
	ModulatorArrayIters mIt;
	BoundElementArrayIters bIt;
} GenesListIterator;

GenesListIterator getGenesListIters(GenesList * g);
void extractInitialQuantities(GenesList * g);

#endif