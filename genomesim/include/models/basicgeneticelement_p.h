//basicgeneticelement_p.h - "Private" structures/functions for basicgeneticelement
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

#ifndef BASIC_GENETIC_ELEMENT_P_H__
#define BASIC_GENETIC_ELEMENT_P_H__
#include "basicgeneticelement.h"
#include "reaction.h"





#define OCC_VEC_PRODUCER_EXISTS_FLAG	1<<0
#define OCC_VEC_ELEM_FLAG(x)	1<<x


typedef struct BitsAndMod {
	occ_bits_rp bits;
	rate_t_rp modifier;
} BitsAndMod;

typedef struct OccupancyVector {
	BitsAndMod *data;
	ulong_type len;
} OccupancyVector;

OccupancyVector emptyOccupancyVector();


struct ProducerElement {
	ElementInfo
	rate_t_rp decayConstant;
	struct {
		rate_t_rp productionConstant;
		pmbPtr produces;
		OccupancyVector occupancies;
		BoundElementArray boundelts;
		ProducerRxns reactions;
	};
	
};

struct ModulatorElement {
	ElementInfo
	rate_t_rp decayConstant;
	BoundElementPtrArray boundelts;
	Reaction * selfDecay;
	ReactionPtrArray bindings;
};

#ifndef DEBUG
inline pmbPtr pmbFromProducer(Producer * p){
	return (pmbPtr) {.producer = p, .species = p->species};
}

inline pmbPtr pmbFromModulator(Modulator * m){
	return (pmbPtr) {.modulator = m, .species = m->species};
}

inline ulong_type getElementID(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->id;
	else if(isModulator(toGet.species)) return toGet.modulator->id;
	else if(isBound(toGet.species)) return toGet.bound->id;
	else return ulong_type_max;
}

inline ulong_type getElementQty(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->qty;
	else if(isModulator(toGet.species)) return toGet.modulator->qty;
	else if(isBound(toGet.species)) return toGet.bound->qty;
	else return ulong_type_max;
}
#else
pmbPtr pmbFromProducer(Producer * p);
pmbPtr pmbFromModulator(Modulator * m);

ulong_type getElementID(pmbPtr toGet);
ulong_type getElementQty(pmbPtr toGet);
#endif

#endif