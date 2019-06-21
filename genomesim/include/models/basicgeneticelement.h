//basicgeneticelement.h - Structure and typedef for genetic elements (producers, modulators, bound)
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
#ifndef BASIC_GENETIC_ELEMENT_H__
#define BASIC_GENETIC_ELEMENT_H__
#include "../globals.h"


struct Reaction;

typedef struct ProducerRxns {
	struct Reaction * production;
	struct Reaction * binding;
	struct Reaction * decay;
	struct Reaction * unbinding;
	struct Reaction * end;
} ProducerRxns;

ProducerRxns nullProducerRxns();
struct ReactionPtrArray;


struct BoundElement;
struct ProducerElement;
struct ModulatorElement;



#define ElementInfo struct { 		\
	char name[MAXGENENAMELENGTH];	\
	species_t species;				\
	ulong_type qty;					\
	ulong_type id;					\
	unsigned char isEnabled;};		



typedef struct pmbPtr{
	union {
		struct ProducerElement * producer;
		struct ModulatorElement * modulator;
		struct BoundElement * bound;
	};
	species_t species;
} pmbPtr;

pmbPtr emptyPmbPtr();



typedef struct BoundElement {
	ElementInfo;
	struct ProducerElement * left;
	struct ModulatorElement * right;
	rate_t_rp assocConstant;
	rate_t_rp dissocConstant;
	effect_t_rp effectStrength;
	struct Reaction * unbinding;
	//relative index in array owned by the producer
	//That is, relativePos = i such that (p.boundElts[i]==this)
	unsigned char producerArrayPos;
	#ifdef BOUNDDECAY
	rate_t_rp decayConstant;
	#endif
} BoundElement;


typedef struct ProducerElement Producer;
typedef struct ModulatorElement Modulator;

typedef  BoundElement * BoundElementPtr;
typedef  Modulator * ModulatorPtr;

typedef struct ProducerElement * ProducerPtr;
struct GenesList;

DECLAREBASICARRAYTYPETYPEDEF(pmbPtr)

DECLAREBASICARRAYTYPETYPEDEF(BoundElement)
DECLAREBASICARRAYTYPETYPEDEF(BoundElementPtr)
DECLAREBASICARRAYTYPETYPEDEF(Producer)
DECLAREBASICARRAYTYPETYPEDEF(ProducerPtr)
DECLAREBASICARRAYTYPETYPEDEF(Modulator)
DECLAREBASICARRAYTYPETYPEDEF(ModulatorPtr)

char isProducer(species_t species);

char isModulator(species_t species);

char isBound(species_t species);

char canDecay(species_t species);

pmbPtr pmbFromBound(BoundElement * b);

ulong_type getProducesID(Producer * p);

void initBasicTF(char *name,Producer * dna, Producer * mess, Modulator * tf,unsigned char enabled);
void initBasicMicro(char *name,Producer *dna,Modulator * miR,unsigned char enabled);	
BoundElement * initConnection(pmbPtr *src, pmbPtr * targ,BoundElement * pos);
void setDefaultProducerQuantity(Producer * p);
void setDefaultModulatorQuantity(Modulator *m);
void setDefaultBoundQuantity(BoundElement * b);
void setProducerQty(Producer * p, unsigned initialQty);
void assignBoundElementsToProducer(Producer *p,const BoundElementArray *bd);
void sortBoundElements(BoundElementArray bd);
void initBoundElementPtrsForModulator(Modulator *m,const BoundElementArray *bd);


ulong_type getIndegree_dna(const Producer *p);
ulong_type getIndegree_mess(const Producer *mess);


#endif
