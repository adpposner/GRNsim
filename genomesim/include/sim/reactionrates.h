// Reactionrates.h - Header for calculating reaction rates and times to events (ttes)
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
#ifndef REACTION_RATES_H__
#define REACTION_RATES_H__
#include "../globals.h"


//Macros for definition of recalc functions 
#define DECLAREEVENTTIMERECALCFUNC(rxn_tp)	 void recalculate ## rxn_tp ## EventTime(struct Reaction * r,unsigned qty); \
 void recalculate ## rxn_tp ## EventTimeSelf(struct Reaction * r,unsigned qty); 

#define DEFINEEVENTTIMERECALCFUNC(rxn_tp)	 void recalculate ## rxn_tp ## EventTime(Reaction * r,unsigned qty){	\
rate_t_rp oldRate = r->reactionRate;	\
rate_t_rp rt = calculate ## rxn_tp ## Rate_(r,qty); \
if (rt<MEPS_RP) setTimeToEvent_(r, INFINITY); \
else if(isinf(r->tte->timeToEvent)) setTimeToEvent_(r, clockTime + getNextExponential( r->reactionRate));\
else{assert( r->tte->timeToEvent >= clockTime);setTimeToEvent_(r, (oldRate/( r->reactionRate))*( r->tte->timeToEvent - clockTime) + clockTime);}} 

#define DEFINEEVENTTIMESELFRECALCFUNC(rxn_tp)  void recalculate ## rxn_tp ## EventTimeSelf(Reaction * r,unsigned qty){	\
rate_t_rp oldRate = r->reactionRate;	\
rate_t_rp rt = r->reactionRate = calculate ## rxn_tp ## Rate_(r, qty); \
if (rt<MEPS_RP) setTimeToEvent_(r, INFINITY); \
else setTimeToEvent_(r,clockTime + getNextExponential(r->reactionRate));}


struct Reaction;
struct ReactionArray;
struct ReactionPtrArray;

struct minHeapNodeArray;

extern time_t_rp clockTime;

void initAllRatesTimes(const struct ReactionArray *arr,struct minHeapNodeArray *mH);
 void updateClock(struct Reaction *r);

DECLAREEVENTTIMERECALCFUNC(Binding)
 void incrementBindingFromLeftEventTime(struct Reaction * r, unsigned qty);

//
 void recalculateProductionEventTimeSelf(struct Reaction *r);
 void recalculateProductionEventTime(struct Reaction *r, unsigned qty);
 void updateProductionEventTime(struct Reaction * r,unsigned molNo,unsigned qty);
DECLAREEVENTTIMERECALCFUNC(Decay)

 void recalculateUnbindingEventTime(struct Reaction *r, unsigned qty);
 //
 void recalculateUnbindingEventTimeSelfNonInf(struct Reaction *r, unsigned qty);
 void setTimeToEventInf(struct Reaction *r);
//
 void setTimeToEventNonZeroRates(struct Reaction * r,rate_t_rp oldRate);
 //
 void setNextEventTimeFromExponential(struct Reaction *r);
 void recalculateEventTime(struct Reaction *r);

#endif