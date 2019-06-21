#include "../include/sim/reactionrates.h"
#include "../include/models/basicgeneticelement.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/randomnumbers.h"
#include "../include/models/reaction.h"
#include "../include/models/minheap.h"
#include "../include/globals.h"


static inline int isSet(const occ_bits_rp_t bits,const unsigned char bitPos){
   
    int elemToCheck = bitPos /64;
    int myInd = bitPos % 64;
    uint64_t lval = bits.x[elemToCheck];
    uint64_t rval = 1 << myInd;
    if (lval & rval) return 1;
    else return 0;
}

static  void calculateDNAProductionModifier_(BoundElementArrayIters bIt, BitsAndMod *bm) {
	assert(bm->bits.x[0] & 1LL);
	numeric_t_rp s = 1.0;
    //bm->modifier = 1.0;
    occ_bits_rp_t tmp = bm->bits;
	ARRAY_TYPE_FOREACH(bIt) {
			if (isSet(tmp,bIt.curr->producerArrayPos))
            s *= bIt.curr->prodEffectStrength;
           
	}
    s = pow(s,TXC_PROD_EXPONENT);
    s = s / (TXC_K_ONE_HALF + s);
    bm->prodMod = s;
}

static  void calculateMessProductionModifier_(BoundElementArrayIters bIt, BitsAndMod *bm) {
	assert(bm->bits.x[0] & 1LL);
	numeric_t_rp pm = 1.0;
    numeric_t_rp dm = 1.0;
    //bm->modifier = 1.0

    occ_bits_rp_t tmp = bm->bits;
	ARRAY_TYPE_FOREACH(bIt) {
		
		if (isSet(tmp,bIt.curr->producerArrayPos)) {
            pm *= bIt.curr->prodEffectStrength;
            assert(bIt.curr->decayEffectStrength > MEPS_RP);
            dm *= bIt.curr->decayEffectStrength;
        }

	}
    pm = powfunc(pm,TXL_PROD_EXPONENT);
    pm = pm / (TXL_K_ONE_HALF + pm);
    dm = powfunc(dm,MESS_DECAY_EXPONENT);
    dm = dm / (MESS_DECAY_K_ONE_HALF + dm);
    bm->prodMod = pm;
    bm->decayMod = dm;
}



// static  rate_t_rp calculateProductionModifiers(BoundElementArrayIters bIt, BitsAndMod bm) {
// 	//If it doesn't exist
// 	assert(bm.bits & 1LL);
// 	rate_t_rp rate = 1.0;
// 	ARRAY_TYPE_FOREACH(bIt) {
// 		bm.bits >>= 1;
// 		if (bm.bits & 1LL) rate *= bIt.curr->effectStrength;
// 	}
// 	return rate;
// }

static  void sumProdRateModifiers_(OccupancyVector * v, unsigned qty) {
	rate_t_rp tmpPM = 0.0;
	BitsAndMod * bm = v->data;
	for (; qty > 0; qty--) {
		tmpPM += (bm++)->prodMod;
	}
	v->prodSum=tmpPM;
}

static  void sumProdDecayRateModifiers_(OccupancyVector * v, unsigned qty) {
	rate_t_rp tmpPM = 0.0;
	rate_t_rp tmpDM = 0.0;
	BitsAndMod * bm = v->data;
	for (; qty > 0; qty--,bm++) {
		tmpPM += bm->prodMod;
		tmpDM += bm->decayMod;
	}
	v->prodSum=tmpPM;
	v->decaySum = tmpDM;
}


static  void calculateDNAProductionRate_(Reaction *p, unsigned qty) {
	Producer * src = p->src;
	p->reactionRate = p->baseRate;
	assert(qty == src->qty);
	BoundElementArrayIters it = getBoundElementArrayIters(&src->boundelts);
	BitsAndMod * oIt = src->occupancies.data;
	unsigned k;
	for (k = qty; k > 0; k--)
		calculateDNAProductionModifier_(it, oIt++);
	sumProdRateModifiers_(&src->occupancies, qty);
	p->reactionRate = p->baseRate * src->occupancies.prodSum;
}


void recalculateMessProductionDecayAfterIncOrDec(Reaction * production, Reaction * decay, unsigned qty){
    	assert(production->src->species == MESSENGER);
    assert(decay->toDecay.species == MESSENGER);
    assert(production->src == decay->toDecay.producer);
    Producer * messenger = production->src;
	assert(qty == messenger->qty);
    sumProdDecayRateModifiers_(&messenger->occupancies,qty);
    production->reactionRate = production->baseRate * messenger->occupancies.prodSum;
    decay->reactionRate = decay->baseRate * messenger->occupancies.decaySum;
}

void calculateMessProductionDecayRate(Reaction *production,Reaction * decay, unsigned qty) {
	assert(production->src->species == MESSENGER);
    assert(decay->toDecay.species == MESSENGER);
    assert(production->src == decay->toDecay.producer);
    Producer * messenger = production->src;
	assert(qty == messenger->qty);
	BoundElementArrayIters it = getBoundElementArrayIters(&messenger->boundelts);
	BitsAndMod * oIt = messenger->occupancies.data;
	unsigned k;
	for (k = qty; k > 0; k--)
		    calculateMessProductionModifier_(it, oIt++);
	sumProdDecayRateModifiers_(&messenger->occupancies,qty);
    production->reactionRate = production->baseRate * messenger->occupancies.prodSum;
    decay->reactionRate = decay->baseRate * messenger->occupancies.decaySum;
}





static  void updateDNAProductionRate_(Reaction *p, int molNo, unsigned qty) {
	Producer *src = p->src;
	p->reactionRate = p->baseRate;
	BitsAndMod * oIt = src->occupancies.data + molNo;
	BoundElementArrayIters it = getBoundElementArrayIters(&src->boundelts);
	calculateDNAProductionModifier_(it, oIt);
	sumProdRateModifiers_(&src->occupancies, qty);
    p->reactionRate = p->baseRate * src->occupancies.prodSum;

}


void updateMessProductionDecayRateAfterProduction(Reaction *production, Reaction * decay,unsigned qty) {
	assert(production->src->species == MESSENGER);
    assert(decay->toDecay.species == MESSENGER);
    assert(production->src == decay->toDecay.producer);
    Producer * messenger = production->src;
	assert(qty == messenger->qty);
	sumProdDecayRateModifiers_(&messenger->occupancies,qty);
	production->reactionRate = production->baseRate * messenger->occupancies.prodSum;
    decay->reactionRate = decay->baseRate * messenger->occupancies.decaySum;

}

void updateMessProductionDecayRateAfterDecay(Reaction *production, Reaction * decay,unsigned qty) {
	assert(production->src->species == MESSENGER);
    assert(decay->toDecay.species == MESSENGER);
    assert(production->src == decay->toDecay.producer);
    Producer * messenger = production->src;
	assert(qty == messenger->qty);
    if(!qty){
	sumProdDecayRateModifiers_(&messenger->occupancies,qty);
	production->reactionRate = production->baseRate * messenger->occupancies.prodSum;
    decay->reactionRate = decay->baseRate * messenger->occupancies.decaySum;}
    else{
        production->reactionRate = 0.0;
        decay->reactionRate = 0.0;
    }
}

void updateMessProductionDecayRate(Reaction *production, Reaction * decay, int molNo, unsigned qty) {
	assert(production->src->species == MESSENGER);
    assert(decay->toDecay.species == MESSENGER);
    assert(production->src == decay->toDecay.producer);
    Producer * messenger = production->src;
	assert(qty == messenger->qty);
	BitsAndMod * oIt = messenger->occupancies.data + molNo;
	BoundElementArrayIters it = getBoundElementArrayIters(&messenger->boundelts);
	calculateDNAProductionModifier_(it, oIt);
	sumProdDecayRateModifiers_(&messenger->occupancies,qty);
	production->reactionRate = production->baseRate * messenger->occupancies.prodSum;
    decay->reactionRate = decay->baseRate * messenger->occupancies.decaySum;
}

static void calculateMessDecayRate_(Reaction * decayRxn,unsigned qty){
    assert(decayRxn->toDecay.species == MESSENGER);
    Reaction * myProduction = decayRxn->toDecay.producer->reactions.production;
    calculateMessProductionDecayRate(myProduction,decayRxn,qty);
}

static void calculateMessProductionRate_(Reaction * prodRxn,unsigned qty){
    assert(prodRxn->src->species == MESSENGER);
    Reaction * myDecay = prodRxn->src->reactions.decay;
    calculateMessProductionDecayRate(prodRxn,myDecay,qty);
}

// static  rate_t_rp calculateMessDecayRate_(Reaction *d, unsigned qty) {
// 	rate_t_rp rate = d->baseRate;
// 	if(!qty)
// 		return d->reactionRate = 0.0;
// 	else if ((d->toDecay.species & MESSENGER)) {
// 		rate = rate * qty;

// 		if(qty != d->toDecay.producer->qty){
// 			printf("%u != %zu\n",qty,d->toDecay.producer->qty);
// 			exit(-203);
// 		}
// 		d->reactionRate = rate; return rate;
// 	} else 
// 	{	
// 		assert(d->toDecay.species & MODULATOR_SPECIES);
// 		rate *= qty;
// 		assert(qty == d->toDecay.modulator->qty);
// 		d->reactionRate = rate; return rate;

// 	}
	
// }

static  void calculateModDecayRate_(Reaction *d, unsigned qty) {
    assert(d->toDecay.species & MODULATOR_SPECIES);
	assert(qty == d->toDecay.modulator->qty);
	rate_t_rp rate = d->baseRate;
	if(!qty)
		d->reactionRate = 0.0;
	 else 
	{	
		rate *= qty;
		d->reactionRate = rate;
	}
	
}



static  void calculateBindingRate_(Reaction *b) {
	rate_t_rp rate = b->baseRate;

	long lq = b->left->qty;
	long rq = b->right->qty;
	long tq = b->target->qty;

	//lq is avail here
	long avail = lq - tq;
	//if(avail<0.0) printf("%f\n",avail);
	if ((avail > 0) && (rq > 0)) {
		rate = rate * (rate_t_rp)(avail * rq);
        b->reactionRate = rate; 
	}
	else {b->reactionRate = 0.0; }
}



static  void calculateUnbindingRate_(Reaction *ub) {
	rate_t_rp rate = ub->baseRate;
	rate *= ub->target->qty;
	ub->reactionRate = rate;
}



static  rate_t_rp calculateReactionRate_(Reaction *r , unsigned qty) {
	unsigned q;
    switch (r->rxn_type) {
	case TFCODING_BINDING:
	case TFNONCODING_BINDING:
	case MESSMIR_BINDING:  calculateBindingRate_(r); return r->reactionRate;
	case TFCODING_UNBINDING:
	case TFNONCODING_UNBINDING:
	case MESSMIR_UNBINDING:  calculateUnbindingRate_(r);return r->reactionRate;
	case MESS_DECAY:  q= r->toDecay.producer->qty;
     calculateMessDecayRate_(r,q);return r->reactionRate;
	case MICRO_DECAY:
    case TF_DECAY: q=r->toDecay.modulator->qty;
     calculateModDecayRate_(r,q);return r->reactionRate;
	case MESS_TRANSCRIPTION:
	case MICRO_TRANSCRIPTION: q=r->src->qty; calculateDNAProductionRate_(r,q);return r->reactionRate;
	case TF_TRANSLATION: q=r->src->qty; calculateMessProductionRate_(r,q);return r->reactionRate;
	default: assert(0);return 0.0;
	}
}

static  void setTimeToEvent_(Reaction *r, time_t_rp t) {
    #ifdef DEBUG
    if(isnan(t)){
        printRxn(stdout,r,1);
        fprintf(stdout,"INVALID TTE NAN\n");
    }
    #endif
	*r->timeToEvent = t;
}

 void setTimeToEventInf(Reaction *r){
	*r->timeToEvent = INFINITY;
}



// void setTimeToEventPossZeroRates(struct Reaction * r);
 void setTimeToEventNonZeroRates(struct Reaction * r, rate_t_rp oldRate) {
	time_t_rp oldTTE = *r->timeToEvent; 
	setTimeToEvent_(r, (oldRate / (r->reactionRate)) * (oldTTE - clockTime) + clockTime);
}

void setTimeToEventPossZeroRates(Reaction * r, rate_t_rp oldRate){
    time_t_rp oldTTE = *r->timeToEvent;
    if(r->reactionRate <= MEPS_RP) setTimeToEventInf(r);
    else if(isinf(oldTTE)) setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
    else {
        setTimeToEventNonZeroRates(r,oldRate);
    }
}

 void setNextEventTimeFromExponential(struct Reaction *r){
	setTimeToEvent_(r, clockTime + getNextExponential(r->reactionRate));
}




 void updateClock(struct Reaction *r) {
	clockTime = *r->timeToEvent;
}





static void calculateEventTimeSelf(Reaction * r, unsigned qty) {
	calculateReactionRate_(r, qty);
	setNextEventTimeFromExponential(r);
}

void initAllRatesTimes(const ReactionArray *arr, minHeap *mH) {
	ReactionArrayIters it = getReactionArrayIters(arr);
	for (it.curr = it.start; it.curr != it.end; it.curr++) {
		it.curr->reactionRate = 0.0;
		calculateEventTimeSelf(it.curr, 0);
	}
	buildMinHeap(mH);
	testMinHeap(*mH);
}

 void updateDNAProductionEventTime(Reaction * r, unsigned qty, unsigned molNo) {
	rate_t_rp oldRate = r->reactionRate;
	/*r->reactionRate = calculate ## rxn_tp ## Rate(r);*/
	updateDNAProductionRate_(r, molNo, qty);
    rate_t_rp newRate = r->reactionRate;
	if (newRate < MEPS_RP) setTimeToEvent_(r, INFINITY);
	else if (oldRate <= MEPS_RP) setTimeToEvent_(r, clockTime + getNextExponential( newRate));
	else {
		assert( *r->timeToEvent >= clockTime); 
		setTimeToEventNonZeroRates(r, oldRate);
	}
}

// void updateMessProductionEventTime(Reaction * r, unsigned qty, unsigned molNo) {
// 	rate_t_rp oldRate = r->reactionRate;
// 	/*r->reactionRate = calculate ## rxn_tp ## Rate(r);*/
// 	updateMessProductionRate_(r, molNo, qty);
//     rate_t_rp newRate = r->reactionRate;
// 	if (newRate < MEPS_RP) setTimeToEvent_(r, INFINITY);
// 	else if (oldRate <= MEPS_RP) setTimeToEvent_(r, clockTime + getNextExponential( newRate));
// 	else {
// 		assert( **r->timeToEvent >= clockTime); 
// 		setTimeToEventNonZeroRates(r, oldRate);
// 	}
// }

 void recalculateBindingEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateBindingRate_(r);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (oldRate <= MEPS_RP) setTimeToEvent_(r, clockTime + getNextExponential(newRate));
	else {
		assert(*r->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r, oldRate);
	}
}




 void recalculateBindingEventTimeSelf(Reaction *r, unsigned qty) {
	calculateBindingRate_(r);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)setTimeToEvent_(r, INFINITY);
	else setTimeToEvent_(r,clockTime + getNextExponential(newRate));
}




 void recalculateUnbindingEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateUnbindingRate_(r);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(*r->timeToEvent)) 
		setNextEventTimeFromExponential(r);
	else {
		assert(*r->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r, oldRate);
	}
}

 void recalculateUnbindingEventTimeSelfNonInf(Reaction *r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateUnbindingRate_(r);
    rate_t_rp newRate = r->reactionRate;
    setTimeToEvent_(r,clockTime + getNextExponential(newRate));
}

//  void recalculateMessDecayEventTime(Reaction * r, unsigned qty) {
// 	rate_t_rp oldRate = r->reactionRate;
// 	calculateMessDecayRate_(r, qty);
//     rate_t_rp newRate = r->reactionRate;
// 	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
// 	else if (isinf(**r->timeToEvent)) setNextEventTimeFromExponential(r);
// 	else {
// 		assert(**r->timeToEvent >= clockTime);
// 		setTimeToEventNonZeroRates(r,oldRate);
// 	}
// }

 void recalculateModDecayEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateModDecayRate_(r, qty);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(*r->timeToEvent)) setNextEventTimeFromExponential(r);
	else {
		assert(*r->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r,oldRate);
	}
}


//  void recalculateMessDecayEventTimeSelf(Reaction *r, unsigned qty) {
// 	rate_t_rp oldRate = r->reactionRate;
// 	calculateMessDecayRate_(r, qty);
//     rate_t_rp newRate = r->reactionRate;
// 	if (newRate <= MEPS_RP)setTimeToEvent_(r, INFINITY);
// 	else setTimeToEvent_(r,clockTime + getNextExponential(newRate));
// }

void recalculateModDecayEventTimeSelf(Reaction *r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateModDecayRate_(r, qty);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)setTimeToEvent_(r, INFINITY);
	else setTimeToEvent_(r,clockTime + getNextExponential(newRate));
}

 void recalculateDNAProductionEventTime(Reaction * r, unsigned qty) {
	rate_t_rp oldRate = r->reactionRate;
	calculateDNAProductionRate_(r, qty);
    rate_t_rp newRate = r->reactionRate;
	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
	else if (isinf(*r->timeToEvent)) setNextEventTimeFromExponential(r);
	else {
		assert(*r->timeToEvent >= clockTime);
		setTimeToEventNonZeroRates(r,oldRate);
	}
}

//  void recalculateMessProductionEventTime(Reaction * r, unsigned qty) {
// 	rate_t_rp oldRate = r->reactionRate;
// 	calculateProductionRate_(r, qty);
//     rate_t_rp newRate = r->reactionRate;
// 	if (newRate <= MEPS_RP)	setTimeToEvent_(r, INFINITY);
// 	else if (isinf(**r->timeToEvent)) setNextEventTimeFromExponential(r);
// 	else {
// 		assert(**r->timeToEvent >= clockTime);
// 		setTimeToEventNonZeroRates(r,oldRate);
// 	}
// }

 void recalculateDNAProductionEventTimeSelf(Reaction *r) {
	setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));
}

//DEFINEEVENTTIMERECALCFUNC(Binding)
//DEFINEEVENTTIMESELFRECALCFUNC(Binding)
//DEFINEEVENTTIMERECALCFUNC(Unbinding)
//DEFINEEVENTTIMESELFRECALCFUNC(Unbinding)
//DEFINEEVENTTIMERECALCFUNC(Decay)
//DEFINEEVENTTIMESELFRECALCFUNC(Decay)
//DEFINEEVENTTIMERECALCFUNC(Production)
//DEFINEEVENTTIMESELFRECALCFUNC(Production)

