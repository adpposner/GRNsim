// reaction.c - Routines for initializing array of reactions and dependency graph
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
#include "../include/models/reaction.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/parameters.h"
#include "../include/models/geneslist.h"
#include "../include/models/geneslist_p.h"
#include "../include/randomnumbers.h"
#include <stdlib.h>
#include "../include/globals.h"
#include <string.h>
#include "../include/models/minheap.h"


void assertRxnProps(Reaction * r,reaction_t rxn_tp,species_t specLeft, species_t specRight){
	assert(r->rxn_type == rxn_tp);
	switch(rxn_tp){
		case PRODUCTION:	if(specLeft) assert(r->src->species & specLeft);
							if(specRight) assert(r->prod.species & specRight);
							break;
		case DECAY:			if(specLeft) assert(r->toDecay.species & specLeft);
						 	break;
		case BINDING:		
		case UNBINDING:		if(specLeft) assert(r->left->species & specLeft);
						 	if(specRight) assert(r->right->species & specRight);
						 	break;
		default: assert(0);
	}
}

// void assertProdByDNA(Reaction *r){assertRxnProps(r, PRODUCTION, DNA, (MESSENGER | MICRO));}
// void assertProdByMess(Reaction *r){ assertRxnProps(r, PRODUCTION, MESSENGER, PROTEIN);}
// void assertDecayByMess(Reaction *r){ assertRxnProps(r, DECAY, MESSENGER, NOSPECIES);}
// void assertDecayByMod(Reaction *r){ assertRxnProps(r, DECAY, MICRO | PROTEIN, NOSPECIES);}
// void assertBindingTFDNA(Reaction *r){assertRxnProps(r, BINDING, DNA, PROTEIN);}
// void assertBindingMessMir(Reaction *r){assertRxnProps(r, BINDING, MESSENGER, MICRO);}
// void assertUnbindingTFDNA(Reaction *r){assertRxnProps(r, UNBINDING, DNA, PROTEIN);}
// void assertUnbindingMessMir(Reaction *r){assertRxnProps(r, UNBINDING, MESSENGER, MICRO);}



 void assertValidBinding(Reaction *r) {
	if (r->left->species & DNA) assert((r->right->species & PROTEIN) && (r->target->species & TFDNA));
	else if (r->left->species & MESSENGER) assert((r->right->species & MICRO) && (r->target->species & MESSMIR));
}

void Reaction_free(Reaction *tofree) {
}

void ReactionPtr_free(ReactionPtr *tofree) {}


DEFINEBASICARRAYTYPEPRIMITIVE(Reaction, Reaction)
DEFINEBASICARRAYTYPEPRIMITIVE(ReactionPtr, ReactionPtr)


const char * event_name(reaction_t e) {
	switch (e) {
	case PRODUCTION:
		return ("PRODUCTION");
	case DECAY:
		return ("DECAY");
	case BINDING:
		return ("BINDING");
	case UNBINDING:
		return ("UNBINDING");
	default:
		return ("");
	}
}

static unsigned currRxnIndex = 0;




//concatenation function for building ReactionArray
ReactionArray reactionArray_cat(ReactionArray * arr, ReactionArray * arr2) {
	ReactionArray toRet;

	toRet = ReactionArray_alloc(arr->length + arr2->length);
	memcpy(toRet.data, arr->data, arr->length * sizeof(*arr->data));
	memcpy(toRet.data + arr->length, arr2->data, arr2->length * (sizeof(*arr2->data)));
	ReactionArray_free(arr);
	ReactionArray_free(arr2);
	return toRet;
}

ReactionPtrArray emptyPtrArray() { return (ReactionPtrArray) {.data = NULL, .length = 0};}

Reaction initProductionReaction(Producer * src) {

	Reaction toRet;
	toRet.id = currRxnIndex++;
	toRet.rxn_type = PRODUCTION;
	toRet.baseRate = src->productionConstant;
	toRet.src = src;
	toRet.prod = src->produces;
	//Don't init dependencies yet
	return toRet;
}

Reaction initBaseDecayReaction(pmbPtr * src) {
	Reaction toRet;
	assert(src->species & (MESSENGER | MICRO | PROTEIN));
	toRet.id = currRxnIndex++;
	toRet.rxn_type = DECAY;
	if (src->species & MESSENGER)
		toRet.baseRate = src->producer->decayConstant;
	else toRet.baseRate = src->modulator->decayConstant;
	toRet.toDecay = *src;
	return toRet;
}

Reaction initBindingReaction(BoundElement * boundElt) {
	Reaction toRet;
	toRet.id  = currRxnIndex++;
	toRet.rxn_type = BINDING;
	toRet.baseRate = boundElt->assocConstant;
	toRet.left = boundElt->left;
	toRet.right = boundElt->right;
	toRet.target = boundElt;
	return toRet;
}

Reaction initUnbindingReaction(BoundElement * boundElt) {
	Reaction toRet;
	toRet.id  = currRxnIndex++;
	toRet.rxn_type = UNBINDING;
	toRet.baseRate = boundElt->dissocConstant;
	toRet.left = boundElt->left;
	toRet.right = boundElt->right;
	toRet.target = boundElt;
	return toRet;
}



static  unsigned numReactionTypeForReactant(const pmbPtr pmb, reaction_t rxn) {
	//Number of reactions is calculated in order to allocate space for the reaction array.
	//Binding/unbinding reactions are lumped with their "left" reactant - i.e. DNA or mRNA
	species_t sp;
	sp = pmb.species;
	if (isEnabled(pmb)) {
		switch (rxn) {
		//# syntheses
		case PRODUCTION: return (isProducer(sp)) ? 1 : 0;
		case DECAY: return (canDecay(sp)) ? 1 : 0;
		case UNBINDING:
		case BINDING: return isBound(sp) ? 1 : 0;
		default: perror("InvalidNumReactionTypeForReactant"); exit(-2313);
		}
	} else return 0;
}

static  unsigned totalReactionsForReactant(const pmbPtr pmb) {
	reaction_t rxntypes[4] = {PRODUCTION, DECAY, BINDING, UNBINDING};
	unsigned toRet;
	toRet = 0;
	int i;
	for (i = 0; i < 4; i++)
		toRet += numReactionTypeForReactant(pmb, rxntypes[i]);
	return toRet;
}

static  ulong_type totalReactionsForGenesList(const GenesList *g) {
	ulong_type totalRxns = 0;
	GenesListIterator gIt = getGenesListIters((GenesList *)g);
	ARRAY_TYPE_FOREACH(gIt.pIt) {
		totalRxns += totalReactionsForReactant(pmbFromProducer(gIt.pIt.curr));
	}
	ARRAY_TYPE_FOREACH(gIt.mIt) {
		totalRxns += totalReactionsForReactant(pmbFromModulator(gIt.mIt.curr));
	}
	ARRAY_TYPE_FOREACH(gIt.bIt) {
		totalRxns += totalReactionsForReactant(pmbFromBound(gIt.bIt.curr));
	}
	return totalRxns;
}

//# of rxns that a producer will hold - this includes all productions, decays,bindings,
//and unbindings which will be init'd in that order
static  ulong_type numRxnsUnderProducer(const pmbPtr pmb) {
	ulong_type toRet = 0;
	assert(isProducer(pmb.species));
	toRet += totalReactionsForReactant(pmb);
	BoundElementArrayIters boundIt = getBoundElementArrayIters(&pmb.producer->boundelts);
	ARRAY_TYPE_FOREACH(boundIt) {
		toRet += totalReactionsForReactant(pmbFromBound(boundIt.curr));
	}
	return toRet;
}

Reaction * initReactionsForBoundElementArray(const BoundElementArray *arr, Reaction * currPos) {
	Reaction * pos = currPos;
	BoundElementArrayIters bIt = getBoundElementArrayIters(arr);
	ARRAY_TYPE_FOREACH(bIt) {
		*pos++ = initBindingReaction(bIt.curr);
	}
	ARRAY_TYPE_FOREACH(bIt) {
		*pos++ = initUnbindingReaction(bIt.curr);
	}
	ulong_type k = 2;
#ifdef BOUNDDECAY
	ARRAY_TYPE_FOREACH(bIt) *pos++ = initBaseDecayReaction(bIt.curr);
	k++;
#endif
	assert(pos == (currPos + (arr->length * k)));
	return pos;
}

Reaction * initReactionsFromProducer(const pmbPtr pmb, ReactionArray * arr, Reaction *currPos) {
	assert(isProducer(pmb.species));
	ulong_type totalRxnsToInit = numRxnsUnderProducer(pmb);
	if (!totalRxnsToInit) return currPos;
	ReactionArrayIters it = getReactionArrayIters(arr);
	assert((currPos - it.start) >= 0);
	assert((it.end - currPos) > 0);
	it.curr = currPos;
	if (numReactionTypeForReactant(pmb, PRODUCTION)) *it.curr++ = initProductionReaction(pmb.producer);
	if (numReactionTypeForReactant(pmb, DECAY)) *it.curr++ = initBaseDecayReaction((pmbPtr *)&pmb);
	it.curr = initReactionsForBoundElementArray(&pmb.producer->boundelts, it.curr);

	assert((it.curr - currPos) == totalRxnsToInit);
	assert(it.curr <= it.end);
	return it.curr;
}

Reaction * initReactionsFromModulator(const pmbPtr pmb, ReactionArray * arr, Reaction * currPos) {
	assert(isModulator(pmb.species));
	ulong_type totalRxnsToInit = numReactionTypeForReactant(pmb, DECAY);
	if(!totalRxnsToInit) return currPos;
	assert(totalRxnsToInit);
	ReactionArrayIters it = getReactionArrayIters(arr);
	assert((currPos - it.start) >= 0);
	assert((it.end - currPos) > 0);
	it.curr = currPos;
	if (numReactionTypeForReactant(pmb, DECAY)) *it.curr++ = initBaseDecayReaction((pmbPtr *)&pmb);
	assert((it.curr - currPos) == totalRxnsToInit);
	assert(it.curr <= it.end);
	return it.curr;
}

int RxnPtrCmpFuncRight(const void *a, const void *b) {
	Reaction ** pa = (Reaction **)a;
	Reaction ** pb = (Reaction **)b;
	reaction_t ra = (*pa)->rxn_type;
	reaction_t rb = (*pb)->rxn_type;
	int dca = 0, dcb = 0;
	if ((ra == PRODUCTION) || (ra == UNBINDING)) dca = 1;
	if ((rb == PRODUCTION) || (rb == UNBINDING)) dcb = 1;
	if ((ra == DECAY) && ((*pa)->toDecay.species == MESSENGER)) dca = 1;
	if ((rb == DECAY) && ((*pb)->toDecay.species == MESSENGER)) dcb = 1;
	if (dca || dcb) return dca - dcb;
	ulong_type ia = (ra == DECAY) ? (*pa)->toDecay.modulator->id :
	                (*pa)->right->id;
	ulong_type ib = (rb == DECAY) ? (*pb)->toDecay.modulator->id :
	                (*pb)->right->id;
	if (ia > ib) return 1;
	else if (ia < ib) return -1;
	else if (ra < rb) return -1;
	else if (ra > rb) return 1;
	else return 0;
}

//Sort on left ID, order Production Binding Decay Unbinding
int RxnCmpFuncLeft(const void *a, const void *b) {
	Reaction * pa = (Reaction *)a;
	Reaction * pb = (Reaction *)b;
	reaction_t ra = pa->rxn_type;
	reaction_t rb = pb->rxn_type;
	int ia = 0, ib = 0, rankA = 0, rankB = 0, dca = 0, dcb = 0;
	switch (ra) {
	case PRODUCTION: ia = pa->src->id; rankA = 0; break;
	case DECAY: if (pa->toDecay.species == MESSENGER) ia = pa->toDecay.producer->id;
		else {ia = pa->toDecay.modulator->id; dca = 1;}
		rankA = 2; break;
	case BINDING: ia = pa->left->id; rankA = 1; break;
	case UNBINDING: ia = pa->left->id; rankA = 3; break;
	default: exit(-134); break;
	}
	switch (rb) {
	case PRODUCTION: ib = pb->src->id; rankB = 0; break;
	case DECAY: if (pb->toDecay.species == MESSENGER) ib = pb->toDecay.producer->id;
		else {ib = pb->toDecay.producer->id; dcb = 1;}
		rankB = 2; break;
	case BINDING: ib = pb->left->id; rankB = 1; break;
	case UNBINDING: ib = pb->left->id; rankB = 3; break;
	default: exit(-252); break;
	}
	if (dca && dcb) return ia - ib;
	else if (dca || dcb) return dca - dcb;
	if (ia > ib) return 1;
	else if (ia < ib) return -1;
	else if (rankA < rankB) return -1;
	else if (rankA > rankB) return 1;
	else return 0;
}

static ReactionArray * reactionArray_init(ulong_type totalLength) {
	ReactionArray * toRet = malloc(sizeof(ReactionArray));
	*toRet = ReactionArray_alloc(totalLength);
	return toRet;
}


static ReactionArray * generateReactionsFromReactants(GenesList *g) {

	ulong_type totalLength = totalReactionsForGenesList(g);
	ReactionArray *	toRet = reactionArray_init(totalLength);
	ReactionArrayIters rIt = getReactionArrayIters(toRet);
	GenesListIterator gIt = getGenesListIters(g);
	rIt.curr = rIt.start;
	ARRAY_TYPE_FOREACH(gIt.pIt) 
		if (isProducerEnabled(gIt.pIt.curr))
			rIt.curr = initReactionsFromProducer((pmbPtr) {.producer = gIt.pIt.curr, .species = gIt.pIt.curr->species},toRet, rIt.curr);
	
	ARRAY_TYPE_FOREACH(gIt.mIt)
		if(isModulatorEnabled(gIt.mIt.curr))
			rIt.curr = initReactionsFromModulator((pmbPtr) {.modulator = gIt.mIt.curr, .species = gIt.mIt.curr->species},toRet, rIt.curr);
	
	assert(rIt.curr == rIt.end);
	qsort(toRet->data, toRet->length, sizeof(*toRet->data), RxnCmpFuncLeft);
	return toRet;
}


static int bindingTargetsAreEqual(const Reaction * a, const  Reaction * b) {
	if (a->left == b->left)
		if (a->right == b->right)
			if (a->target == b->target)
				return 1;
	return 0;
}

static void getUnbindingForBinding(Reaction * binding, const Reaction * unbindingStart, const Reaction * unbindingEnd) {
	assert(binding->rxn_type == BINDING);
	if (unbindingStart == unbindingEnd)
		return;
	const Reaction * unbindingCurr;
	for (unbindingCurr = unbindingStart; unbindingCurr != unbindingEnd; unbindingCurr++) {
		assert(unbindingCurr->rxn_type == UNBINDING);
		if (bindingTargetsAreEqual(binding, unbindingCurr))
		{binding->target->unbinding = (Reaction *)unbindingCurr; return;}
	}
	assert(0);
}

static void initUnbindingsForProducer(Producer *p) {
	const Reaction *ubStart = p->reactions.unbinding;
	const Reaction *ubEnd = p->reactions.end;
	ReactionArrayIters rIt = {.start = p->reactions.binding, .end = p->reactions.decay, .curr = p->reactions.binding};
	ARRAY_TYPE_FOREACH(rIt) getUnbindingForBinding(rIt.curr, ubStart, ubEnd);
}



ReactionPtrArray getReactionPtrsSortedRight(ReactionArray *r) {
	ReactionPtrArray toRet = ReactionPtrArray_alloc(r->length);
	ReactionArrayIters rIt = getReactionArrayIters(r);
	ReactionPtrArrayIters pIt = getReactionPtrArrayIters(&toRet);
	pIt.curr = pIt.start;
	ARRAY_TYPE_FOREACH(rIt) {
		*pIt.curr++ = rIt.curr;
	}
	assert(pIt.curr == pIt.end);
	qsort(toRet.data, toRet.length, sizeof(*toRet.data), RxnPtrCmpFuncRight);
	return toRet;
}

static void ensureModulatorsAreSorted(ModulatorArray *m) {
	ModulatorArrayIters it = getModulatorArrayIters(m);
	ulong_type lastID = 0;
	ARRAY_TYPE_FOREACH(it) {
		assert(it.curr->id > lastID);
		lastID = it.curr->id;
	}
}
static void ensureProducersAreSorted(ProducerArray *m) {
	ProducerArrayIters it = getProducerArrayIters(m);
	ulong_type lastID = it.start->id;
	it.curr = it.start + 1;
	for (; it.curr != it.end; it.curr++) {
		assert(it.curr->id > lastID);
		lastID = it.curr->id;
	}
}



static void initDependenciesForModulators(GenesList *g, ReactionArray *rxns) {
	ReactionPtrArray arr = getReactionPtrsSortedRight(rxns);
	ensureModulatorsAreSorted(&g->modulators);
	ModulatorArrayIters mIt = getModulatorArrayIters(&g->modulators);
	ReactionPtrArrayIters rIt = getReactionPtrArrayIters(&arr);
	rIt.curr = rIt.start;
	ReactionPtr * startPos, * endPos;
	ulong_type nBindingReactions = 0;
	ARRAY_TYPE_FOREACH(mIt) {
		//loop once for count
		nBindingReactions = 0;
		assert((*rIt.curr)->rxn_type == DECAY);
		mIt.curr->selfDecay = *rIt.curr++;
		startPos = rIt.curr;
		while ((*rIt.curr)->rxn_type == BINDING) {
			nBindingReactions++;
			rIt.curr++;
		}
		endPos = rIt.curr;
		mIt.curr->bindings = ReactionPtrArray_alloc(nBindingReactions);
		ReactionPtrArrayIters rpIt = getReactionPtrArrayIters(&mIt.curr->bindings);
		rpIt.curr = rpIt.start;
		for (rIt.curr = startPos; rIt.curr != endPos; rIt.curr++) {
			*rpIt.curr++ = *rIt.curr;
			assert((*rIt.curr)->rxn_type == BINDING);
		}
		assert(rpIt.curr == rpIt.end);
	}
	ReactionPtrArray_free(&arr);
}


static void initDependenciesForProducers(GenesList *g, ReactionArray *rxns) {
	ensureProducersAreSorted(&g->producers);
	ProducerArrayIters pIt = getProducerArrayIters(&g->producers);
	ReactionArrayIters rIt = getReactionArrayIters(rxns);
	ARRAY_TYPE_FOREACH(pIt) {
		pIt.curr->reactions = nullProducerRxns();
		assert(rIt.curr->rxn_type == PRODUCTION);
		pIt.curr->reactions.production = rIt.curr++;
		if (rIt.curr->rxn_type == BINDING) {
			pIt.curr->reactions.binding = rIt.curr++;
			assert(rIt.curr->left == pIt.curr);
		}
		while (rIt.curr->rxn_type == BINDING) {
			assert(rIt.curr->left == pIt.curr);
			rIt.curr++;
		}
		if (pIt.curr->species == MESSENGER) {
			assert(rIt.curr->left == pIt.curr);
			assert(rIt.curr->rxn_type == DECAY);
			pIt.curr->reactions.decay = rIt.curr++;

		}
		if (rIt.curr->rxn_type == UNBINDING) {
			assert(rIt.curr->left == pIt.curr);
			pIt.curr->reactions.unbinding = rIt.curr++;
		}
		while (rIt.curr->rxn_type == UNBINDING) {
			assert(rIt.curr->left == pIt.curr);
			rIt.curr++;
		}
		pIt.curr->reactions.end = rIt.curr;
		if (!pIt.curr->reactions.unbinding) pIt.curr->reactions.unbinding = pIt.curr->reactions.end;
		if (!pIt.curr->reactions.decay) pIt.curr->reactions.decay = pIt.curr->reactions.unbinding;
		if (!pIt.curr->reactions.binding) pIt.curr->reactions.binding = pIt.curr->reactions.decay;
		initUnbindingsForProducer(pIt.curr);
	}
}

static void generateReactionDependencies(GenesList *g, ReactionArray *rxns) {
	initDependenciesForProducers(g, rxns);
	initDependenciesForModulators(g, rxns);
}

void generateReactionsAndDependencies(GenesList *g, ReactionArray **rxns) {
	*rxns = generateReactionsFromReactants(g);
	generateReactionDependencies(g, *rxns);
}



////END REACTION GENERATION ROUTINES ////////

///BEGIN UTILITIES	///////////////////////////
void printRxn(FILE * fh, const Reaction * rxn, const int wDeps) {
#define PRINTRXNBUFSIZE	1000
	char buf[PRINTRXNBUFSIZE] = {0};
	memset(buf, 0, PRINTRXNBUFSIZE);
	int cx = 0;
	int n = 0;



	n += snprintf(buf + cx, PRINTRXNBUFSIZE - cx, "ID = %u\t", rxn->id); if (n >= 0)cx += n;

	switch (rxn->rxn_type) {
	case PRODUCTION: cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s: ", "PRODUCTION");

		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "SRCTYPE = %s SRC ID = %zu, PROD ID = %zu SRCQTY= %zu, PRODQTY = %zu", species_names(rxn->src->species),
		               rxn->src->id, getElementID(rxn->prod), rxn->src->qty, getElementQty(rxn->prod));
		break;
	case DECAY:	cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s: ", "DECAY");
		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "SPTYPE = %s ID = %zu QTY = %zu ",
		               species_names(rxn->toDecay.species), getElementID(rxn->toDecay), getElementQty(rxn->toDecay));
		break;
	case UNBINDING: cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s", "UN");
	case BINDING: cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s: ", "BINDING");

		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "LEFT: %zu, RIGHT:%s, %zu, TARG: %zu QTY = %zu ", rxn->left->id, (rxn->right->species == PROTEIN) ? "PROT" : "MIR", rxn->right->id, rxn->target->id, rxn->target->qty);
		break;

	default: fprintf(fh, "UNKNOWN REACTION TYPE %s:%d TYPE = %d", __FILE__, __LINE__, rxn->rxn_type); //exit(-5);
	}
	cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, " TTE: %f RATE: %f", rxn->tte->timeToEvent, rxn->reactionRate);
	cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, " DEPS: ");

	fprintf(fh, "\t%s\n", buf);


}

void printAllRxns(FILE * fh, ReactionArray *r) {
	ReactionArrayIters it = getReactionArrayIters(r);
	ARRAY_TYPE_FOREACH(it) printRxn(fh, it.curr, 0);
}

