#ifndef REACTION_RATES_H__
#define REACTION_RATES_H__
#include "../globals.h"
struct Reaction;
struct ReactionArray;
struct ReactionPtrArray;

struct minHeapSoA;

extern time_t_rp clockTime;

void initAllRatesTimes(const struct ReactionArray *arr,struct minHeapSoA *mH);
 void updateClock(struct Reaction *r);

// #define DECLAREEVENTTIMERECALCFUNC(rxn_tp)	 void recalculate ## rxn_tp ## EventTime(struct Reaction * r,unsigned qty); \
//  void recalculate ## rxn_tp ## EventTimeSelf(struct Reaction * r,unsigned qty); 

// #define DEFINEEVENTTIMERECALCFUNC(rxn_tp)	 void recalculate ## rxn_tp ## EventTime(Reaction * r,unsigned qty){	\
// rate_t_rp oldRate = r->reactionRate;	\
// /*r->reactionRate = calculate ## rxn_tp ## Rate(r);*/ \
// rate_t_rp rt = calculate ## rxn_tp ## Rate_(r,qty); \
// if (rt<MEPS_RP) setTimeToEvent_(r, INFINITY); \
// else if(isinf(r->tte->timeToEvent)) setTimeToEvent_(r, clockTime + getNextExponential( r->reactionRate));\
// else{assert( r->tte->timeToEvent >= clockTime);setTimeToEvent_(r, (oldRate/( r->reactionRate))*( r->tte->timeToEvent - clockTime) + clockTime);}} 

// #define DEFINEEVENTTIMESELFRECALCFUNC(rxn_tp)  void recalculate ## rxn_tp ## EventTimeSelf(Reaction * r,unsigned qty){	\
// rate_t_rp oldRate = r->reactionRate;	\
// rate_t_rp rt = r->reactionRate = calculate ## rxn_tp ## Rate_(r, qty); \
// if (rt<MEPS_RP) setTimeToEvent_(r, INFINITY); \
// else setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));}

//DECLAREEVENTTIMERECALCFUNC(Binding)
void recalculateBindingEventTime(struct Reaction * r, unsigned qty);
void recalculateBindingEventTimeSelf(struct Reaction * r, unsigned qty);
 void incrementBindingFromLeftEventTime(struct Reaction * r, unsigned qty);
//DECLAREEVENTTIMERECALCFUNC(Production)
 void recalculateDNAProductionEventTimeSelf(struct Reaction *r);
 void recalculateDNAProductionEventTime(struct Reaction *r, unsigned qty);
 void updateDNAProductionEventTime(struct Reaction * r,unsigned molNo,unsigned qty);

 //void recalculateMessDecayEventTime(struct Reaction * r, unsigned qty);

 void recalculateModDecayEventTime(struct Reaction * r, unsigned qty);

// void recalculateMessDecayEventTimeSelf(struct Reaction *r, unsigned qty);
void updateMessProductionDecayRateAfterDecay(struct Reaction *production,struct  Reaction * decay,unsigned qty);
void recalculateModDecayEventTimeSelf(struct Reaction *r, unsigned qty);
void updateMessProductionDecayRateAfterProduction(struct Reaction *production,struct Reaction * decay,unsigned qty);


//DECLAREEVENTTIMERECALCFUNC(Unbinding)
 void recalculateUnbindingEventTime(struct Reaction *r, unsigned qty);
 void recalculateUnbindingEventTimeSelfNonInf(struct Reaction *r, unsigned qty);
 void setTimeToEventInf(struct Reaction *r);

 void setTimeToEventPossZeroRates(struct Reaction * r,rate_t_rp oldRate);
 void setTimeToEventNonZeroRates(struct Reaction * r,rate_t_rp oldRate);
 void setNextEventTimeFromExponential(struct Reaction *r);

void recalculateMessProductionDecayAfterIncOrDec(struct Reaction * production,struct Reaction * decay, unsigned qty);
void calculateMessProductionDecayRate(struct Reaction *production,struct Reaction * decay, unsigned qty);
void updateMessProductionDecayRate(struct Reaction *production,struct  Reaction * decay, int molNo, unsigned qty);
#endif