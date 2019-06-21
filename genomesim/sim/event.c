// Event.c - contains all minheap update funcs and clock updates
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
#include "../include/sim/event.h"
#include "../include/globals.h"
#include "../include/randomnumbers.h"
#include "../include/sim/reactionrates.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/models/reaction.h"
#include <float.h>
#ifdef DEBUG
#include <math.h>
#endif
//#include <mathimf.h>
#include "../include/models/minheap.h"


synthTotals totalSynths = {0, 0, 0};

ulong_type totalProductions = 0;
ulong_type totalDecays = 0;
ulong_type totalBindings = 0;
ulong_type totalUnbindings = 0;
ulong_type totalForcedUnbindings = 0;
ulong_type totalForcedDecays = 0;

//These allow us to record events in a similar manner for both
//miRNA present and absent situations by only tracking things that are
//non-miRNA-related
long commonBirthDeaths = 0;
long commonBirthDeathsLast = 0;


//Experimental macros to make recalc & update functions - not currently used
/*
#define DEFINERECALCUPDATEFUNC(rxn_tp) static  void recalculate ## rxn_tp ##AndPushDown(Reaction * r,minHeap *mH,unsigned qty) {\
rate_t_rp oldRate = r->reactionRate; recalculate ## rxn_tp ##EventTime(r,qty); \
assert(r->reactionRate <= oldRate); pushDown(mH,r->tte->index);} \
static  void recalculate ## rxn_tp ##AndPushUp(Reaction *r,minHeap *mH,unsigned qty){\
rate_t_rp oldRate = r->reactionRate; recalculate ## rxn_tp ##EventTime(r,qty); \
assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);}

#define DEFINERECALCSELFUPDATEFUNC(rxn_tp) static  void recalculateSelf ## rxn_tp ##AndUpdateHeap(Reaction * r,minHeap *mH,unsigned qty) {\
recalculate ## rxn_tp ##EventTimeSelf(r,qty); updateRoot(mH);}
*/

//
static  void recalculateBindingAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->tte->index);
}

//
static  void incrementBindingFromLeftAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void recalculateBindingAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void recalculateSelfBindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateBindingEventTimeSelf(r, qty); updateRoot(mH);
}

//
static  void recalculateProductionAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateProductionEventTime(r, qty);
	//Depending on the exponent used for miRNAs, small variations can occur
	//here, causing the need for fabs
	assert(fabs(r->reactionRate - oldRate) < 0.01f); pushDown(mH, r->tte->index);
}

//
static  void incrementProductionAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	updateProductionEventTime(r, qty, qty - 1);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void recalculateSelfProductionAndUpdateHeap(Reaction *r, minHeap *mH) {
	recalculateProductionEventTimeSelf(r); updateRoot(mH);
}



//
static  void recalculateDecayAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateDecayEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->tte->index);
}

//
static  void recalculateModDecayAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateDecayEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void incrementMessDecayAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateDecayEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void recalculateSelfDecayAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateDecayEventTimeSelf(r, qty); updateRoot(mH);
}

//
static  void recalculateUnbindingAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateUnbindingEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->tte->index);
}

//
static  void recalculateUnbindingAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateUnbindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->tte->index);
}

//
static  void recalculateSelfUnbindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	if (qty)
	{recalculateUnbindingEventTimeSelfNonInf(r, qty); updateRoot(mH);}
	else {
		setTimeToEventInf(r); pushDown(mH, r->tte->index);
	}
}




//
static  void updateProductionAndUpdateHeap(Reaction * production, minHeap *mH, unsigned qty, unsigned molNo) {
	rate_t_rp oldRate = production->reactionRate;
	assert(qty == production->src->qty);
	updateProductionEventTime(production, qty, molNo);
	if (oldRate < production->reactionRate)
		pushUp(mH, production->tte->index);
	else
		pushDown(mH, production->tte->index);
}

//
static  void recalculateRightAndPushUp(Modulator * right, minHeap *mH, const unsigned qty) {
	recalculateModDecayAndPushUp(right->selfDecay, mH, right->qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushUp(*it.curr, mH, qty);
	}
}

//
static  void recalculateRightAndPushDown(Modulator * right, minHeap *mH, unsigned qty) {
	recalculateDecayAndPushDown(right->selfDecay, mH, qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushDown(*it.curr, mH, qty);
	}
}

//
static  void recalculateRightAndUpdateHeap(Modulator * right, minHeap *mH, unsigned qty) {
	recalculateDecayEventTime(right->selfDecay, qty);
	update(mH, right->selfDecay->tte->index);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingEventTime(*it.curr, qty);
		update(mH, (*it.curr)->tte->index);
	}
}



//following mess synth
//
static  void incrementLeftAndPushUp(ProducerRxns prodrxns, minHeap *mH, unsigned qty) {

	incrementProductionAndPushUp(prodrxns.production, mH, qty);
	Reaction * curr;
	for (curr = prodrxns.binding; curr != prodrxns.decay; curr++)
		incrementBindingFromLeftAndPushUp(curr, mH, qty);
	incrementMessDecayAndPushUp(prodrxns.decay, mH, qty);

}

//Following mess decay
//
static  void recalculateLeftAndPushDown(ProducerRxns prodrxns, minHeap *mH, unsigned qty) {

	recalculateProductionAndPushDown(prodrxns.production, mH, qty);
	Reaction * curr;
	for (curr = prodrxns.binding; curr != prodrxns.decay; curr++)
		recalculateBindingAndPushDown(curr, mH, qty);
	recalculateDecayAndPushDown(prodrxns.decay, mH, qty);
	for (curr = prodrxns.unbinding; curr != prodrxns.end; curr++)
		recalculateUnbindingAndPushDown(curr, mH, qty);

}



//
static  void decrementProducerQty(Producer *p) {
	assert(p->qty > 0);
	p->qty--;
}



//
unsigned setFirstUnoccupiedRandomStart(Producer * p,
                                       const unsigned char bitPos) {
	const occ_bits_rp bits = 1LL << bitPos;
	unsigned leftQty = p->qty;
	assert(leftQty <= MAX_MRNA_QTY);
	unsigned offset = getNextBindingPosition(leftQty);

	int i;
	for (i = 0; i < leftQty; i++) {
		if (p->occupancies.data[(i + offset) % leftQty].bits & bits)
			continue;
		p->occupancies.data[(i + offset) % leftQty].bits |= bits;
		return (i + offset) % leftQty;
	}
	perror("can't find bits to set"); exit(-8421);
}

//
unsigned clearFirstOccupiedRandomStart(Producer * p,
                                       const unsigned char bitPos) {
	const occ_bits_rp bits = 1LL << bitPos;
	unsigned leftQty = p->qty;
	assert(leftQty <= MAX_MRNA_QTY);
	unsigned offset = getNextBindingPosition(leftQty);
	int i;
	//keep moving to next molecule until bitCondition is not equal to condition on molecule
	for (i = 0; i < leftQty; i++) {
		if (!(p->occupancies.data[(i + offset) % leftQty].bits & bits))
			continue;
		p->occupancies.data[(i + offset) % leftQty].bits &= ~bits;
		return (i + offset) % leftQty;
	}
	perror("can't find bits to clear"); exit(-8422);
}

//
unsigned executeBindingReaction(Reaction * binding) {
	binding->target->qty++; binding->right->qty--;
	totalBindings++;

	if (binding->left->species & CODING) {
		commonBirthDeaths++;
	}

	return setFirstUnoccupiedRandomStart(binding->left, binding->target->producerArrayPos);
}

//
unsigned executeUnbindingReaction(Reaction * unbinding) {
	unbinding->target->qty--; unbinding->right->qty++;
	totalUnbindings++;

	if (unbinding->left->species & CODING) {
		commonBirthDeaths++;
	}

	return clearFirstOccupiedRandomStart(unbinding->left, unbinding->target->producerArrayPos);
	
}



//
static  void executeForcedUnbindingReaction(BoundElement * pb, minHeap *mH) {
	pb->qty--; pb->right->qty++;
	recalculateUnbindingAndPushDown(pb->unbinding, mH, pb->qty);
	recalculateRightAndUpdateHeap(pb->right, mH, pb->right->qty);
	totalForcedUnbindings++;
}

#include <math.h>

//
void executeForcedDecayReaction(BoundElement * pb, minHeap *mH) {
	pb->qty--;
	recalculateUnbindingAndPushDown(pb->unbinding, mH, pb->qty);
	totalForcedDecays++;
}



//
void executeModDecayReaction(Reaction * decay) {
	decay->toDecay.modulator->qty--;
	totalDecays++;
	if (decay->toDecay.species & (PROTEIN | MESSENGER))
		commonBirthDeaths++;
}

//
static  void swapOccupancies(BitsAndMod * curr, BitsAndMod *last) {
	BitsAndMod tmp;
	assert(curr->bits & 1);
	assert(last->bits & 1);
	if (curr == last) return;
	tmp = *curr;
	*curr = *last;
	*last = tmp;
}

//
rate_t_rp executeMessDecayReaction(Reaction * decay, minHeap *mH) {

	Producer * mess = decay->toDecay.producer;
	//Choose a number corresponding to the index to erase
	unsigned ind = getNextBindingPosition(mess->qty);
	decrementProducerQty(mess);
	//Move the one to be destroyed to the end
	BitsAndMod *last = mess->occupancies.data + mess->qty;
	BitsAndMod *curr = mess->occupancies.data + ind;
	swapOccupancies(curr, last);
	totalDecays++;
	BoundElementArrayIters bIt = getBoundElementArrayIters(&mess->boundelts);
	occ_bits_rp lastBits = last->bits;
	assert(lastBits & 1);
	ARRAY_TYPE_FOREACH(bIt) {
		lastBits >>= 1;
		if (lastBits & 1) {
#ifdef MESSENGERDECAYFORCESUNBINDING
#ifndef MESSENGERDECAYFORCESMICRODECAY
			executeForcedUnbindingReaction(bIt.curr, mH);
#else
#error "cannot set both MESSENGERDECAYFORCESUNBINDING and MESSENGERDECAYFORCESMICRODECAY"
#endif
#elif defined MESSENGERDECAYFORCESMICRODECAY
			executeForcedDecayReaction(bIt.curr, mH);
#else
#error "must define either MESSENGERDECAYFORCESUNBINDING or MESSENGERDECAYFORCESMICRODECAY"
#endif
		}
	}
	rate_t_rp toRet = last->modifier;
	last->modifier = 0.0;
	commonBirthDeaths++;
	return toRet;

}



//
void executeModProductionReaction(Reaction * pr) {

	pr->prod.modulator->qty++;
	totalProductions++;
	if (pr->prod.species & (PROTEIN | MESSENGER))
		commonBirthDeaths++;
}

//
rate_t_rp executeMessProductionReaction(Reaction * pr) {
	Producer * p = pr->prod.producer;
	p->occupancies.data[p->qty] = (BitsAndMod) {.bits = OCC_VEC_PRODUCER_EXISTS_FLAG, .modifier = 1.0};
	p->qty++;
	totalProductions++;
	commonBirthDeaths++;
	return 1.0;
}

//
void executeReaction(minHeap *mH) {
	Reaction *toExec = mH->data[0].rxn;
	updateClock(toExec);
	unsigned molNo;
	reaction_t rxn_type = toExec->rxn_type;
	if (rxn_type == BINDING) {
		molNo = executeBindingReaction(toExec);
		recalculateSelfBindingAndUpdateHeap(toExec, mH, 0);
		updateProductionAndUpdateHeap(toExec->left->reactions.production, mH, toExec->left->qty, molNo);
		recalculateUnbindingAndPushUp(toExec->target->unbinding, mH, toExec->target->qty);
		recalculateRightAndPushDown(toExec->right, mH, toExec->right->qty);
	}
	else if (rxn_type == UNBINDING) {
		molNo = executeUnbindingReaction(toExec);
		recalculateSelfUnbindingAndUpdateHeap(toExec, mH, toExec->target->qty);
		updateProductionAndUpdateHeap(toExec->left->reactions.production, mH, toExec->left->qty, molNo);
		recalculateRightAndPushUp(toExec->right, mH, toExec->right->qty);
	}
	else if (rxn_type == PRODUCTION) {
		if (toExec->prod.species == MESSENGER) {
			executeMessProductionReaction(toExec);
			totalSynths.nMessTranscription++;
			recalculateSelfProductionAndUpdateHeap(toExec, mH);
			incrementLeftAndPushUp(toExec->prod.producer->reactions, mH, toExec->prod.producer->qty);
		} else {
			executeModProductionReaction(toExec);
			if (toExec->src->species == MESSENGER)
				totalSynths.nTranslation++;
			else
				totalSynths.nMicroTranscription++;
			recalculateSelfProductionAndUpdateHeap(toExec, mH);
			recalculateRightAndPushUp(toExec->prod.modulator, mH, toExec->prod.modulator->qty);
		}
	} else if (rxn_type == DECAY) {

		if (toExec->toDecay.species == MESSENGER) {
			executeMessDecayReaction(toExec, mH);
			recalculateSelfDecayAndUpdateHeap(toExec, mH, toExec->toDecay.producer->qty);
			recalculateLeftAndPushDown(toExec->toDecay.producer->reactions, mH, toExec->toDecay.producer->qty);
		} else {
			executeModDecayReaction(toExec);
			recalculateSelfDecayAndUpdateHeap(toExec, mH, toExec->toDecay.modulator->qty);
			recalculateRightAndPushDown(toExec->toDecay.modulator, mH, toExec->toDecay.modulator->qty);
		}
	}


#ifdef TESTASSERTIONS
if (testMinHeap(*mH)) {
	printRxn(stdout, toExec, 1);
	exit(95);
}
#endif

}

void printTotals() {
	printf("prods="ltf_rp"\td="ltf_rp"\tb="ltf_rp"\tu="ltf_rp", fu="ltf_rp" fd="ltf_rp"\n", totalProductions, totalDecays, totalBindings, 
	totalUnbindings, totalForcedUnbindings, totalForcedDecays);
	totalProductions = totalDecays = totalBindings = totalUnbindings = totalForcedUnbindings = totalForcedDecays = commonBirthDeaths = 0;
}



