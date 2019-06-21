#ifndef REACTION_H__
#define REACTION_H__

#include "../globals.h"
#include "basicgeneticelement.h"
#include "geneslist.h"

typedef enum reaction_t {EMPTYRXN=0,
MESS_TRANSCRIPTION = 1<<0, MICRO_TRANSCRIPTION = 1<<1, TF_TRANSLATION = 1<<2,
ANY_TRANSCRIPTION = (1<<0) | (1<<1), ANY_PRODUCTION = (1<<0) | (1<<1) | (1<<2),
MESS_DECAY=1<<3, MICRO_DECAY = 1<<4, TF_DECAY = 1<<5, MOD_DECAY = (1<<4) | (1<<5),
ANY_DECAY = (1<<3) | (1<<4) | (1<<5),
MESSMIR_BINDING=1<<6, TFCODING_BINDING = 1<<7, TFNONCODING_BINDING = 1<< 8,
TFDNA_BINDING = (1<<7) | (1<<8), ANY_BINDING = (1<<6) | (1<<7) | (1<<8),
MESSMIR_UNBINDING=1<<9, TFCODING_UNBINDING = 1<<10, TFNONCODING_UNBINDING = 1<<11,
TFDNA_UNBINDING = (1<<10) | (1<<11), ANY_UNBINDING = (1<<9) | (1<<10) | (1<<11)} 
reaction_t;

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
	double * timeToEvent;
    ulong_type mHIndex;
	rate_t_rp reactionRate;
}	Reaction;

typedef struct Reaction * ReactionPtr;
DECLAREBASICARRAYTYPEPRIMITIVE(ReactionPtr, ReactionPtr)


DECLAREBASICARRAYTYPEPRIMITIVE(Reaction, Reaction)

//ReactionArray generateReactionsFromReactants(GenesList *g);
//void generateReactionDependencies(GenesList *g, ReactionArray *rxns);
//void initDependenciesForModulators(GenesList *g, ReactionArray *rxns);
//void initDependenciesForProducers(GenesList *g, ReactionArray *rxns);
void generateReactionsAndDependencies(GenesList *g, ReactionArray **rxns);

//After a TF is produced, we simply keep the rate the same so we get a new exponential
void assertRxnProps(Reaction * r,reaction_t rxn_tp,species_t specLeft, species_t specRight);



void assertValidBinding(Reaction *r);

 ReactionArrayIters getDependenciesForDecayOfMessenger(Reaction *r);
 ReactionArrayIters getDependenciesForProductionOfMessenger(Reaction *r);


void printRxn(FILE * fh,const Reaction * rxn,const int wDeps);
void printAllRxns(FILE * fh,ReactionArray *r);
void writeReactions_JSON(ReactionArray arr,const char * filename);


#endif