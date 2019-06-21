// Reactionrates.c -Calculating reaction rates and times to events (ttes)
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
#include "../include/sim/reactionrates.h"
#include "../include/models/basicgeneticelement.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/randomnumbers.h"
#include "../include/models/reaction.h"
#include "../include/models/minheap.h"
#include <math.h>
#include <float.h>
#include "../include/globals.h"




//
static  void calculateProductionModifier(BoundElementArrayIters bIt, BitsAndMod *bm) {
	assert(bm->bits & 1LL);
	bm->modifier = 1.0;
	occ_bits_rp checkbit = 1 << 0;
	ARRAY_TYPE_FOREACH(bIt) {
		checkbit <<= 1;
		if (bm->bits & checkbit) bm->modifier *= bIt.curr->effectStrength;
	}

}



//
static  rate_t_rp sumRateModifiers(OccupancyVector * v, unsigned qty) {
	rate_t_rp toRet = 0.0;
	BitsAndMod * bm = v->data;
	for (; qty > 0; qty--) {
		toRet += (bm++)->modifier;
	}
	return toRet;
}

//
static  rate_t_rp calculateProductionRate_(Reaction *p, unsigned qty) {
	Producer * src = p->src;
	p->reactionRate = p->baseRate;
	assert(qty == src->qty);
	BoundElementArrayIters it = getBoundElementArrayIters(&src->boundelts);
	BitsAndMod * oIt = src->occupancies.data;
	unsigned k;
	for (k = qty; k > 0; k--)
		calculateProductionModifier(it, oIt++);
	p->reactionRate *=  sumRateModifiers(&src->occupancies, qty);
	return p->reactionRate;
}

//
static  rate_t_rp updateProductionRate_(Reaction *p, int molNo, unsigned qty) {
	Producer *src = p->src;
	p->reactionRate = p->baseRate;
	BitsAndMod * oIt = src->occupancies.data + molNo;
	BoundElementArrayIters it = getBoundElementArrayIters(&src->boundelts);
	calculateProductionModifier(it, oIt);
	p->reactionRate *= sumRateModifiers(&src->occupancies, qty);
	return p->reactionRate;
}

//
static  rate_t_rp calculateDecayRate_(Reaction *d, unsigned qty) {
	rate_t_rp rate = d->baseRate;
	if(!qty)
		return d->reactionRate = 0.0;
	else if ((d->toDecay.species & MESSENGER)) {
		rate = rate * qty;

		if(qty != d->toDecay.producer->qty){
			printf("%u != %zu\n",qty,d->toDecay.producer->qty);
			exit(-203);
		}
		d->reactionRate = rate; return rate;
	} else 
	{	
		assert(d->toDecay.species & MODULATOR_SPECIES);
		rate *= qty;
		assert(qty == d->toDecay.modulator->qty);
		d->reactionRate = rate; return rate;

	}
	
}



//
static inline rate_t_rp calculateBindingRate_(Reaction *b, unsigned qty) {
	rate_t_rp rate = b->baseRate;

	int  avail = (b->left->qty - b->target->qty)*b->right->qty;
	if (avail > 0) {
		 rate = (rate_t_rp)avail * rate; b->reactionRate = rate; return rate;

	}
	else {b->reactionRate = 0.0; return b->reactionRate;}
}



//
static inline  rate_t_rp calculateUnbindingRate_(Reaction *ub, unsigned qty) {
	rate_t_rp rate = ub->baseRate;
	rate *= ub->target->qty;
	ub->reactionRate = rate;
	return rate;
}

//
static  rate_t_rp calculateReactionRate_(Reaction *r , unsigned qty) {
	switch (r->rxn_type) {
	case BINDING: return r->reactionRate = calculateBindingRate_(r, qty);
	case UNBINDING: return r->reactionRate = calculateUnbindingRate_(r, qty);
	case PRODUCTION: return r->reactionRate = calculateProductionRate_(r, r->src->qty);
	case DECAY: return r->reactionRate = calculateDecayRate_(r, (r->toDecay.species == MESSENGER) ?
		r->toDecay.producer->qty : r->toDecay.modulator->qty);
	default: return 0.0;
	}
}

//
static  void setTimeToEvent_(Reaction *r, time_t_rp t) {
	r->tte->timeToEvent = t;
}

//
void setTimeToEventInf(Reaction *r){
	r->tte->timeToEvent = INFINITY;
}

//
void setTimeToEventNonZeroRates(struct Reaction * r, rate_t_rp oldRate) {
	time_t_rp oldTTE = r->tte->timeToEvent; 
	setTimeToEvent_(r, (oldRate / (r->reactionRate)) * (oldTTE - clockTime) + clockTime);
}

//
void setNextEventTimeFromExponential(struct Reaction *r){
	setTimeToEvent_(r, clockTime + getNextExponential(r->reactionRate));
}

//
 void updateClock(struct Reaction *r) {
	clockTime = r->tte->timeToEvent;
}

//
static void calculateEventTimeSelf(Reaction * r, unsigned qty) {
	calculateReactionRate_(r, qty);
	setNextEventTimeFromExponential(r);
}

//
void initAllRatesTimes(const ReactionArray *arr, minHeap *mH) {
	ReactionArrayIters it = getReactionArrayIters(arr);
	for (it.curr = it.start; it.curr != it.end; it.curr++) {
		it.curr->reactionRate = 0.0;
		calculateEventTimeSelf(it.curr, 0);
	}
	buildMinHeap(mH);
	testMinHeap(*mH);
}

//
 void updateProductionEventTime(Reaction * r, unsigned qty, unsigned molNo) {
	rate_t_rp oldRate = r->reactionRate;
	/*r->reactionRate = calculate ## rxn_tp ## Rate(r);*/
	rate_t_rp rt = updateProductionRate_(r, molNo, qty);
	if (rt < MEPS_RP) setTimeToEvent_(r, INFINITY);
	else if (oldRate <= MEPS_RP) setTimeToEvent_(r, clockTime + getNextExponential( r->reactionRate));
	else {
		assert( r->tte->timeToEvent >= clockTime); 
		setTimeToEventNonZeroRates(r, oldRate);	}
}

//
 void recalculateBindingEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateBindingRate_(r, qty);
	if (rt <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(r->tte->timeToEvent)) setNextEventTimeFromExponential(r);
	else {
		assert(r->tte->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r, oldRate);
	}
}

//
 void recalculateBindingEventTimeSelf(Reaction *r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateBindingRate_(r, qty);
	if (rt <= MEPS_RP)setTimeToEvent_(r, INFINITY);
	else setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
}

//
 void recalculateUnbindingEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateUnbindingRate_(r, qty);
	if (rt <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(r->tte->timeToEvent)) 
		setNextEventTimeFromExponential(r);
	else {
		assert(r->tte->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r, oldRate);
	}
}

//
 void recalculateUnbindingEventTimeSelfNonInf(Reaction *r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateUnbindingRate_(r, qty);
	setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
}

//
 void recalculateDecayEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateDecayRate_(r, qty);
	if (rt <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(r->tte->timeToEvent)) setNextEventTimeFromExponential(r);
	else {
		assert(r->tte->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r,oldRate);
	}
}

//
 void recalculateDecayEventTimeSelf(Reaction *r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateDecayRate_(r, qty);
	if (rt <= MEPS_RP)setTimeToEvent_(r, INFINITY);
	else setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
}

//
 void recalculateProductionEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	rate_t_rp rt = calculateProductionRate_(r, qty);
	if (rt <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(r->tte->timeToEvent)) setNextEventTimeFromExponential(r);
	else {
		assert(r->tte->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r,oldRate);
	}
}

//
 void recalculateProductionEventTimeSelf(Reaction *r) {
	setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
}

