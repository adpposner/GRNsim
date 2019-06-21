#include "../include/sim/event.h"
#include "../include/globals.h"
#include "../include/randomnumbers.h"
#include "../include/sim/reactionrates.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/models/reaction.h"
#include <float.h>
#include <mathimf.h>
#include "../include/models/minheap.h"


synthTotals totalSynths = {0};


//These allow us to record events in a similar manner for both
//miRNA present and absent situations by only tracking things that are
//non-miRNA-related
long commonBirthDeaths = 0;
long commonBirthDeathsLast = 0;

// #define DEFINERECALCUPDATEFUNC(rxn_tp) static  void recalculate ## rxn_tp ##AndPushDown(Reaction * r,minHeap *mH,unsigned qty) {\
// rate_t_rp oldRate = r->reactionRate; recalculate ## rxn_tp ##EventTime(r,qty); \
// assert(r->reactionRate <= oldRate); pushDown(mH,r->mH);} \
// static  void recalculate ## rxn_tp ##AndPushUp(Reaction *r,minHeap *mH,unsigned qty){\
// rate_t_rp oldRate = r->reactionRate; recalculate ## rxn_tp ##EventTime(r,qty); \
// assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);}

// #define DEFINERECALCSELFUPDATEFUNC(rxn_tp) static  void recalculateSelf ## rxn_tp ##AndUpdateHeap(Reaction * r,minHeap *mH,unsigned qty) {\
// recalculate ## rxn_tp ##EventTimeSelf(r,qty); updateRoot(mH);}

static  void recalculateBindingAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->mHIndex);
}

static  void incrementBindingFromLeftAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
}

static  void recalculateBindingAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateBindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
}



static  void recalculateSelfBindingAndUpdateHeap_(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateBindingEventTimeSelf(r, qty); updateRoot(mH);
}

static  void recalculateSelfMessMirBindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateSelfBindingAndUpdateHeap_(r,mH,qty);
}

static  void recalculateSelfTFDNABindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateSelfBindingAndUpdateHeap_(r,mH,qty);
}


// static  void incrementMessProductionDecayAndPushUp(ProducerRxns * prodRxns, minHeap *mH, unsigned int qty) {
// 	rate_t_rp oldProdRate = prodRxns->production->reactionRate;
//     rate_t_rp oldDecayRate = prodRxns->decay->reactionRate;
//     recalculateMessProductionDecayAfterIncOrDec(prodRxns->production,prodRxns->decay,qty);
//     if(!qty){
//         setTimeToEventInf(prodRxns->production);
//         setTimeToEventInf(prodRxns->decay);
//     }else{
	
//     }
//     assert(0);
// 	//assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
// }

static  void recalculateSelfDNAProductionAndUpdateHeap(Reaction *r, minHeap *mH) {
    recalculateDNAProductionEventTimeSelf(r); updateRoot(mH);
}



static  void recalculateModDecayAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateModDecayEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->mHIndex);
}

static  void recalculateModDecayAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateModDecayEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
}

static  void recalculateMessDecayAndPushUp(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateModDecayEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
}




static  void recalculateSelfModDecayAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateModDecayEventTimeSelf(r, qty); updateRoot(mH);
}


static  void recalculateUnbindingAndPushDown(Reaction *r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateUnbindingEventTime(r, qty);
	assert(r->reactionRate <= oldRate); pushDown(mH, r->mHIndex);
}

static  void recalculateUnbindingAndPushUp_(Reaction * r, minHeap *mH, unsigned int qty) {
	rate_t_rp oldRate = r->reactionRate;
	recalculateUnbindingEventTime(r, qty);
	assert(r->reactionRate >= oldRate); pushUp(mH, r->mHIndex);
}

static  void recalculateMessMirUnbindingAndPushUp(Reaction * r, minHeap *mH, unsigned int qty){
	recalculateUnbindingAndPushUp_(r,mH,qty);
}
static  void recalculateTFDNAUnbindingAndPushUp(Reaction * r, minHeap *mH, unsigned int qty){
	recalculateUnbindingAndPushUp_(r,mH,qty);
}

static  void recalculateSelfUnbindingAndUpdateHeap_(Reaction *r, minHeap *mH, unsigned qty) {
	if (qty)
	{recalculateUnbindingEventTimeSelfNonInf(r, qty); updateRoot(mH);}
	else {
		setTimeToEventInf(r); pushDown(mH, r->mHIndex);
	}
}

static  void recalculateSelfMessMirUnbindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateSelfUnbindingAndUpdateHeap_(r,mH,qty);
}

static  void recalculateSelfTFDNAUnbindingAndUpdateHeap(Reaction *r, minHeap *mH, unsigned qty) {
	recalculateSelfUnbindingAndUpdateHeap_(r,mH,qty);
}



// //If the executed reaction is a production, then the rate is unchanged,
// //so we don't need to recalculate it
// static  void recalculateProductionSelfAndUpdateHeap(Reaction * production,minHeap *mH,unsigned qty){
// 	recalculateProductionEventTimeSelf(production,qty);
// 	assert(mH->data->rxn == production);
// 	updateRoot(mH);
// }

// static  void recalculateProductionAndUpdateHeap(Reaction * production, minHeap *mH, unsigned qty) {
// 	rate_t_rp oldRate = production->reactionRate;
// 	recalculateProductionEventTime(production, qty);
// 	if (oldRate < production->reactionRate)
// 		pushUp(mH, production->mHIndex);
// 	else
// 		pushDown(mH, production->mHIndex);
// }

static  void updateDNAProductionAndUpdateHeap_(Reaction * production, minHeap *mH, unsigned qty, unsigned molNo) {
	rate_t_rp oldRate = production->reactionRate;
	assert(qty == production->src->qty);
	updateDNAProductionEventTime(production, qty, molNo);
	if (oldRate < production->reactionRate)
		pushUp(mH, production->mHIndex);
	else
		pushDown(mH, production->mHIndex);
}


static void recalculateSelfTFProductionAndUpdateHeap(Reaction * r, minHeap * mH){
    setNextEventTimeFromExponential(r);
    updateRoot(mH);
}

static  void updateDNAProductionAndUpdateHeap(Reaction * production, minHeap *mH, unsigned qty, unsigned molNo){
    updateDNAProductionAndUpdateHeap_(production,mH,qty,molNo);
}



static  void recalculateModAfterUnbindingAndPushUp(Modulator * right, minHeap *mH, const unsigned qty) {
	recalculateModDecayAndPushUp(right->selfDecay, mH, right->qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushUp(*it.curr, mH, qty);
	}
}

static  void recalculateModAfterProductionAndPushUp(Modulator * right, minHeap *mH, const unsigned qty) {
	recalculateModDecayAndPushUp(right->selfDecay, mH, right->qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushUp(*it.curr, mH, qty);
	}
}

static  void recalculateModAfterBindingAndPushDown(Modulator * right, minHeap *mH, unsigned qty) {
	recalculateModDecayAndPushDown(right->selfDecay, mH, qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushDown(*it.curr, mH, qty);
	}
}

static  void recalculateModAfterDecayAndPushDown(Modulator * right, minHeap *mH, unsigned qty) {
	recalculateModDecayAndPushDown(right->selfDecay, mH, qty);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingAndPushDown(*it.curr, mH, qty);
	}
}


static  void recalculateRightAndUpdateHeap(Modulator * right, minHeap *mH, unsigned qty) {
	recalculateModDecayEventTime(right->selfDecay, qty);
	update(mH, right->selfDecay->mHIndex);
	ReactionPtrArrayIters it = getReactionPtrArrayIters(&right->bindings);
	ARRAY_TYPE_FOREACH(it) {
		recalculateBindingEventTime(*it.curr, qty);
		update(mH, (*it.curr)->mHIndex);
	}
}

// static  void recalculateProductionRateAndPushUp(Producer *p,minHeap *mH,unsigned qty){
// 	recalculateProduction
// }





//following mess synth
static  void incrementMessProductionDecayAndPushUp(ProducerRxns prodrxns, minHeap *mH, unsigned qty) {
	//incrementProductionDecayAndPushUp(prodrxns, mH, qty);
	rate_t_rp oldDecayRate = prodrxns.decay->reactionRate;
    rate_t_rp oldProdRate = prodrxns.production->reactionRate;
    updateMessProductionDecayRateAfterProduction(prodrxns.production,prodrxns.decay,qty);
    setTimeToEventPossZeroRates(prodrxns.production,oldProdRate);
    pushUp(mH,prodrxns.production->mHIndex);
    setTimeToEventPossZeroRates(prodrxns.decay,oldDecayRate);
    pushUp(mH,prodrxns.decay->mHIndex);

	// for(curr=prodrxns.unbinding;curr!=prodrxns.end;curr++)
	// 	recalculateUnbindingAndUpdateHeap(curr, mH);

}

static void decrementMessProductionDecayAndPushDown(ProducerRxns prodrxns, minHeap * mH, unsigned qty){
    rate_t_rp oldProdRate = prodrxns.production->reactionRate;
    updateMessProductionDecayRateAfterDecay(prodrxns.production,prodrxns.decay,qty);
    if(!qty){
        setTimeToEventInf(prodrxns.production);
        pushDown(mH,prodrxns.production->mHIndex);
        setTimeToEventInf(prodrxns.decay);
        pushDown(mH,prodrxns.decay->mHIndex);
    }else{
        setNextEventTimeFromExponential(prodrxns.decay);
        updateRoot(mH);
        setTimeToEventPossZeroRates(prodrxns.production,oldProdRate);
        pushDown(mH,prodrxns.production->mHIndex);
    }
}



//Following mess decay
static  void recalculateMessAfterDecayAndPushDown(ProducerRxns prodrxns, minHeap *mH, unsigned qty) {
    decrementMessProductionDecayAndPushDown(prodrxns,mH,qty);
    //recalculateProductionDecayAndPushDown(prodrxns.production, mH, qty);
	Reaction * curr;
	for (curr = prodrxns.binding; curr != prodrxns.decay; curr++)
		recalculateBindingAndPushDown(curr, mH, qty);
//	recalculateMessDecayAndPushDown(prodrxns.decay, mH, qty);
	for (curr = prodrxns.unbinding; curr != prodrxns.end; curr++)
		recalculateUnbindingAndPushDown(curr, mH, qty);

}

static void recalculateMessAfterProductionAndPushUp(ProducerRxns prodrxns, minHeap *mH, unsigned qty){
    incrementMessProductionDecayAndPushUp(prodrxns,mH,qty);
        Reaction * curr;
	for (curr = prodrxns.binding; curr != prodrxns.decay; curr++)
		incrementBindingFromLeftAndPushUp(curr, mH, qty);
}


static void updateMessProductionDecayAndUpdateAfterBinding(ProducerRxns prodrxns, minHeap * mH, unsigned qty,unsigned molNo){
    //production always goes down, decay always goes up
    Reaction * production = prodrxns.production;
    rate_t_rp oldProdRate = production->reactionRate;
    Reaction * decay = prodrxns.decay;
    rate_t_rp oldDecayRate = decay->reactionRate;
    updateMessProductionDecayRate(production,decay,molNo,qty);
    setTimeToEventPossZeroRates(production,oldProdRate);
    pushUp(mH,production->mHIndex);
    pushDown(mH,production->mHIndex);
    setTimeToEventPossZeroRates(decay,oldDecayRate);
    pushUp(mH,decay->mHIndex);
}

static void updateMessProductionDecayAndUpdateAfterUnbinding(ProducerRxns prodrxns, minHeap * mH, unsigned qty,unsigned molNo){
    //production always goes down, decay always goes up
    Reaction * production = prodrxns.production;
    rate_t_rp oldProdRate = production->reactionRate;
    Reaction * decay = prodrxns.decay;
    rate_t_rp oldDecayRate = decay->reactionRate;
    updateMessProductionDecayRate(production,decay,molNo,qty);
    setTimeToEventPossZeroRates(production,oldProdRate);
    pushUp(mH,production->mHIndex);
    setTimeToEventPossZeroRates(decay,oldDecayRate);
    pushUp(mH,decay->mHIndex);
    pushDown(mH,decay->mHIndex);
}


// #define OCCUPANCY_ALTERING_FUNCS
// #ifdef	OCCUPANCY_ALTERING_FUNCS
static  void incrementProducerQty(Producer * p) {
	//because of array zero-indexing we want to set the occupancy before
	//incrementing qty
	p->qty++;
}

static  void decrementProducerQty(Producer *p) {
	//must do something before
	assert(p->qty > 0);
	p->qty--;
}

static inline int isBitSet(Producer * p,const int molNo,const unsigned char bitPos){
    assert(bitPos != 0);
    int elemToCheck = bitPos /64;
    int myInd = bitPos % 64;
    uint64_t lval = p->occupancies.data[molNo].bits.x[elemToCheck];
    uint64_t rval = 1 << myInd;
    if (lval & rval) return 1;
    else return 0;
}


static inline int isBitSetAtPos(const occ_bits_rp_t bits,const unsigned char bitPos){
   
    int elemToCheck = bitPos /64;
    int myInd = bitPos % 64;
    uint64_t lval = bits.x[elemToCheck];
    uint64_t rval = 1 << myInd;
    if (lval & rval) return 1;
    else return 0;
}



static inline void toggleBitAt(Producer * p,const int molNo,const unsigned char bitPos){
    
    int elemToCheck = bitPos /64;
    int myInd = bitPos % 64;
    
    uint64_t rval = 1 << myInd;
    
    p->occupancies.data[molNo].bits.x[elemToCheck] ^= rval;
}




// #endif
unsigned setFirstUnoccupiedRandomStart(Producer * p,
                                       const unsigned char bitPos) {
	const occ_bits_rp_t bits = {.x = {1 << bitPos, 0, 0 ,0}};
	unsigned leftQty = p->qty;
	assert(leftQty <= MAX_MRNA_QTY);
	unsigned offset = getNextBindingPosition(leftQty);

	int i;
    int molNo;
	for (i = 0; i < leftQty; i++) {
        molNo = (i+offset) % leftQty;
        if(isBitSet(p,molNo,bitPos)) continue;
        else{
            toggleBitAt(p,molNo,bitPos);
    		return molNo;
        }
	}
	perror("can't find bits to set"); exit(-8421);
}

unsigned clearFirstOccupiedRandomStart(Producer * p,
                                       const unsigned char bitPos) {
	const occ_bits_rp_t bits = {.x = {1 << bitPos, 0, 0 ,0}};
	unsigned leftQty = p->qty;
	assert(leftQty <= MAX_MRNA_QTY);
	unsigned offset = getNextBindingPosition(leftQty);
	int i;
    int molNo;
	//keep moving to next molecule until bitCondition is not equal to condition on molecule
	for (i = 0; i < leftQty; i++) {
        molNo = (i+offset) % leftQty;
		if (!isBitSet(p,molNo,bitPos))
			continue;
		toggleBitAt(p,molNo,bitPos);
		return molNo;
	}
	perror("can't find bits to clear"); exit(-8422);
}


// unsigned executeBindingReaction(Reaction * binding) {
// 	binding->target->qty++; binding->right->qty--;
// 	totalSynths.totalb

// 	if (binding->left->species & CODING) {
// 		commonBirthDeaths++;
// 	}

// 	return setFirstUnoccupiedRandomStart(binding->left, binding->target->producerArrayPos);
// 	//binding->left->occupancies.data[pos].modifier *= binding->target->effectStrength;
// }

unsigned executeMessMirBindingReaction(Reaction * binding) {
		binding->target->qty++; binding->right->qty--;
		totalSynths.totalMessMirBindings++;

	return setFirstUnoccupiedRandomStart(binding->left, binding->target->producerArrayPos);
}
unsigned executeTFDNABindingReaction(Reaction * binding) {
		binding->target->qty++; binding->right->qty--;
		totalSynths.totalTFDNABindings++;

	return setFirstUnoccupiedRandomStart(binding->left, binding->target->producerArrayPos);
}

unsigned executeUnbindingReaction_(Reaction * unbinding){
	unbinding->target->qty--; unbinding->right->qty++;
	return clearFirstOccupiedRandomStart(unbinding->left, unbinding->target->producerArrayPos);
}

unsigned executeTFDNAUnbindingReaction(Reaction * unbinding) {

	totalSynths.totalTFDNAUnBindings++;
	return executeUnbindingReaction_(unbinding);
}

unsigned executeMessMirUnbindingReaction(Reaction * unbinding) {
	totalSynths.totalMessMirUnBindings++;
	return executeUnbindingReaction_(unbinding);

}




static  void executeForcedUnbindingReaction(BoundElement * pb, minHeap *mH) {
	pb->qty--; pb->right->qty++;
	recalculateUnbindingAndPushDown(pb->unbinding, mH, pb->qty);
	recalculateRightAndUpdateHeap(pb->right, mH, pb->right->qty);
	totalSynths.totalForcedMirUnbindings++;
}

#include <math.h>

void setZeroRate(Reaction * r, minHeap *mH) {
	r->reactionRate = 0.0;
	*r->timeToEvent = INFINITY;
}

// void decrementRateAndUpdate(Reaction * ub, minHeap * mH) {
// 	rate_t_rp oldRate = ub->reactionRate;
// 	ub->reactionRate -= ub->baseRate;
// 	setTimeToEventNonZeroRates(ub, oldRate);
// }

void executeForcedDecayReaction(BoundElement * pb, minHeap *mH) {
	pb->qty--;
	recalculateUnbindingAndPushDown(pb->unbinding, mH, pb->qty);
	totalSynths.totalForcedMirDecays++;
}



void executeModDecayReaction(Reaction * decay) {
	//Modulator * de = decay->toDecay.modulator;
	decay->toDecay.modulator->qty--;
	totalSynths.totalModDecays++;
}

static  void swapOccupancies(BitsAndMod * curr, BitsAndMod *last) {
	BitsAndMod tmp;
	assert(curr->bits.x[0] & 1);
	assert(last->bits.x[0] & 1);
	if (curr == last) return;
	tmp = *curr;
	*curr = *last;
	*last = tmp;
	//can change rate here
}

unsigned selectMessToDestroy(Producer * messr){
    rate_t_rp dS = messr->occupancies.decaySum;
    effect_t_rp sel;
    randomUniformNumberArray(1,0.0,dS,&sel);
    BitsAndMod * bm = messr->occupancies.data;
    unsigned ind;
    for(ind=0;ind<messr->qty;ind++,bm++){
        sel -= bm->decayMod;
        if (sel <= MEPS_RP)
            return ind;
    }
    fprintf(stderr,"UNABLE TO FIND MESS TO DESTROY, sel at end == %f\n",sel);
    exit(-23424);
}

void executeMessDecayReaction(Reaction * decay, minHeap *mH) {
totalSynths.totalMessDecays++;
	Producer * mess = decay->toDecay.producer;
	//Choose a number corresponding to the index to erase
	unsigned ind = selectMessToDestroy(mess);
	// int lastInd = mess->qty-1;
	decrementProducerQty(mess);
	//Move the one to be destroyed to the end
	BitsAndMod *last = mess->occupancies.data + mess->qty;
	BitsAndMod *curr = mess->occupancies.data + ind;
	swapOccupancies(curr, last);
	
	BoundElementArrayIters bIt = getBoundElementArrayIters(&mess->boundelts);
	occ_bits_rp_t lastBits = last->bits;
	assert(isBitSetAtPos(lastBits,0));
    int i=1;
	ARRAY_TYPE_FOREACH(bIt) {
		if(isBitSetAtPos(lastBits,bIt.curr->producerArrayPos)){
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
	last->bits.x[0] = last->bits.x[1]  = last->bits.x[2] = last->bits.x[3] =0 ;
	last->decayMod=0.0;
	last->prodMod=0.0;
}




void executeModProductionReaction(Reaction * pr) {
	//Modulator * m = pr->prod.modulator;
	pr->prod.modulator->qty++;
	if (pr->prod.species == PROTEIN)
		totalSynths.totalTranslations++;
	else
		totalSynths.totalNonCodingTranscription++;
}

rate_t_rp executeMessProductionReaction(Reaction * pr) {
	Producer * p = pr->prod.producer;
	p->occupancies.data[p->qty] = (BitsAndMod) {.bits = {.x = {1<<0,0,0,0}}, .prodMod = 1.0,.decayMod=1.0};
	p->qty++;
	totalSynths.totalCodingTranscription++;
	return 1.0;
}

void executeDecayReaction(Reaction * pr) {}

void executeReaction(minHeap *mH) {
	Reaction *toExec = mH->rxnArr[0];
	updateClock(toExec);
	//printf("%f\n",clockTime);
	//printRxn(stdout, toExec, 0);
	unsigned molNo;
	reaction_t rxn_type = toExec->rxn_type;
	//break down all types
	if (rxn_type & MESSMIR_BINDING) {
		totalSynths.totalMessMirBindings++;
		molNo = executeMessMirBindingReaction(toExec);
		recalculateSelfMessMirBindingAndUpdateHeap(toExec, mH, 0);
		updateMessProductionDecayAndUpdateAfterBinding(toExec->left->reactions, mH, toExec->left->qty, molNo);
		recalculateMessMirUnbindingAndPushUp(toExec->target->unbinding, mH, toExec->target->qty);
		recalculateModAfterBindingAndPushDown(toExec->right, mH, toExec->right->qty);
	} else if (rxn_type & MESSMIR_UNBINDING){
		totalSynths.totalMessMirUnBindings++;
		molNo = executeMessMirUnbindingReaction(toExec);
		recalculateSelfMessMirUnbindingAndUpdateHeap(toExec, mH, toExec->target->qty);
		updateMessProductionDecayAndUpdateAfterUnbinding(toExec->left->reactions, mH, toExec->left->qty, molNo);
		recalculateModAfterUnbindingAndPushUp(toExec->right, mH, toExec->right->qty);
	} else if (rxn_type & TFDNA_BINDING){
		totalSynths.totalTFDNABindings++;
		molNo = executeTFDNABindingReaction(toExec);
		recalculateSelfTFDNABindingAndUpdateHeap(toExec, mH, 0);
		updateDNAProductionAndUpdateHeap(toExec->left->reactions.production, mH, toExec->left->qty, molNo);
		recalculateTFDNAUnbindingAndPushUp(toExec->target->unbinding, mH, toExec->target->qty);
		recalculateModAfterBindingAndPushDown(toExec->right, mH, toExec->right->qty);
	}else if (rxn_type & TFDNA_UNBINDING){
		totalSynths.totalTFDNAUnBindings++;
		molNo = executeTFDNAUnbindingReaction(toExec);
		recalculateSelfTFDNAUnbindingAndUpdateHeap(toExec, mH, toExec->target->qty);
		updateDNAProductionAndUpdateHeap(toExec->left->reactions.production, mH, toExec->left->qty, molNo);
		recalculateModAfterUnbindingAndPushUp(toExec->right, mH, toExec->right->qty);
	}else if (rxn_type & MESS_DECAY){
		totalSynths.totalMessDecays++;
			executeMessDecayReaction(toExec, mH);
			recalculateMessAfterDecayAndPushDown(toExec->toDecay.producer->reactions, mH, toExec->toDecay.producer->qty);
	}else if (rxn_type & MOD_DECAY){
		totalSynths.totalModDecays++;
			executeModDecayReaction(toExec);
			recalculateSelfModDecayAndUpdateHeap(toExec, mH, toExec->toDecay.modulator->qty);
			recalculateModAfterDecayAndPushDown(toExec->toDecay.modulator, mH, toExec->toDecay.modulator->qty);
	}else if (rxn_type & MESS_TRANSCRIPTION){
		totalSynths.totalCodingTranscription++;
		executeMessProductionReaction(toExec);
		recalculateSelfDNAProductionAndUpdateHeap(toExec, mH);
		recalculateMessAfterProductionAndPushUp(toExec->prod.producer->reactions, mH, toExec->prod.producer->qty);
	}else if (rxn_type & MICRO_TRANSCRIPTION){
			executeModProductionReaction(toExec);
			totalSynths.totalNonCodingTranscription++;
			recalculateSelfDNAProductionAndUpdateHeap(toExec, mH);
			recalculateModAfterProductionAndPushUp(toExec->prod.modulator, mH, toExec->prod.modulator->qty);

	}else if (rxn_type & TF_TRANSLATION){
			executeModProductionReaction(toExec);
			totalSynths.totalTranslations++;
			recalculateSelfTFProductionAndUpdateHeap(toExec, mH);
			recalculateModAfterProductionAndPushUp(toExec->prod.modulator, mH, toExec->prod.modulator->qty);

	}else{assert(0);}

//recalculateEventTimesWithDependencies(mH,toExec);
//updateDependencies(mH, toExec);


//printRxn(stdout, toExec, 1);
//updateDependencies(mH, toExec);


#ifdef TESTASSERTIONS
if (testMinHeap(*mH)) {
	printRxn(stdout, toExec, 1);
//		generate_from_mh(*mH, toExec);
	exit(95);
}
#endif

}

void printTotals() {
	fprintf(stdout," totalCodingTranscription = "ltf_rp" totalNonCodingTranscription = "ltf_rp
	 " totalTranslations = "ltf_rp" totalMessDecays = "ltf_rp
	 " totalForcedMirDecays = "ltf_rp" totalForcedMirUnbindings = "ltf_rp" totalModDecays = "ltf_rp
	 " totalMessMirBindings = "ltf_rp" totalTFDNABindings = "ltf_rp" totalMessMirUnBindings = "ltf_rp
	 " totalTFDNAUnBindings = "ltf_rp"\n",	
	 totalSynths.totalCodingTranscription,totalSynths.totalNonCodingTranscription,
	totalSynths.totalTranslations,totalSynths.totalMessDecays,
	totalSynths.totalForcedMirDecays,totalSynths.totalForcedMirUnbindings,
	totalSynths.totalModDecays,totalSynths.totalMessMirBindings,
	totalSynths.totalTFDNABindings,totalSynths.totalMessMirUnBindings,
	totalSynths.totalTFDNAUnBindings);
	
 totalSynths.totalCodingTranscription = totalSynths.totalNonCodingTranscription = totalSynths.totalTranslations = 
 totalSynths.totalMessDecays = totalSynths.totalForcedMirDecays = totalSynths.totalForcedMirUnbindings = 
 totalSynths.totalModDecays = totalSynths.totalMessMirBindings = totalSynths.totalTFDNABindings = 
 totalSynths.totalMessMirUnBindings = totalSynths.totalTFDNAUnBindings = 0;
}



