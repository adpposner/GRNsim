//reaction.h - Structures for reaction object
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
#ifndef REACTION_H__
#define REACTION_H__

#include "../globals.h"
#include "basicgeneticelement.h"
#include "geneslist.h"

typedef enum reaction_t {EMPTYRXN=0,PRODUCTION =1,DECAY=2,BINDING=4,UNBINDING=8} reaction_t;

typedef struct productionReaction {
	Producer *src;
	pmbPtr prod;
} productionReaction;

typedef struct decayReaction {
	pmbPtr toDecay;
} decayReaction;

typedef struct bindingReaction {
	Producer * left;
	Modulator * right;
	BoundElement * target;
} bindingReaction;

typedef union reactionComponents {
	productionReaction pr;
	decayReaction de;
	bindingReaction bu;
} reactionComponents;


struct minHeapNode;

typedef struct Reaction {
	unsigned id;
	reaction_t rxn_type;
	rate_t_rp baseRate;
	union {
		struct { 
			Producer *src;
			pmbPtr prod;
		};
		struct {
			pmbPtr toDecay;
		};
		struct {Producer *left;
			Modulator *right;
			BoundElement *target;
		};
	};
	struct minHeapNode *tte;
	rate_t_rp reactionRate;
}	Reaction;

typedef struct Reaction * ReactionPtr;
DECLAREBASICARRAYTYPEPRIMITIVE(ReactionPtr, ReactionPtr)


DECLAREBASICARRAYTYPEPRIMITIVE(Reaction, Reaction)


void generateReactionsAndDependencies(GenesList *g, ReactionArray **rxns);


void assertRxnProps(Reaction * r,reaction_t rxn_tp,species_t specLeft, species_t specRight);


void assertValidBinding(Reaction *r);

 ReactionArrayIters getDependenciesForDecayOfMessenger(Reaction *r);
 ReactionArrayIters getDependenciesForProductionOfMessenger(Reaction *r);


void printRxn(FILE * fh,const Reaction * rxn,const int wDeps);
void printAllRxns(FILE * fh,ReactionArray *r);
void writeReactions_JSON(ReactionArray arr,const char * filename);


#endif