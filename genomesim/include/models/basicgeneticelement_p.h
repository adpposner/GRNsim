#ifndef BASIC_GENETIC_ELEMENT_P_H__
#define BASIC_GENETIC_ELEMENT_P_H__
#include "basicgeneticelement.h"
#include "reaction.h"





// #define OCC_VEC_PRODUCER_EXISTS_FLAG	1<<0
// #define OCC_VEC_ELEM_FLAG(x)	1<<x


typedef struct BitsAndMod {
	occ_bits_rp_t bits;
	rate_t_rp prodMod;
	rate_t_rp decayMod;
} BitsAndMod;

typedef struct OccupancyVector {
	BitsAndMod *data;
	rate_t_rp prodSum;
	rate_t_rp decaySum;
	ulong_type len;
} OccupancyVector;

OccupancyVector emptyOccupancyVector();

// typedef struct OccupancyVectorIters {
// 	occ_bits_rp * const occStart,* const occEnd;
// 	rate_t_rp * const modStart,* const modEnd;
// 	occ_bits_rp * occCurr;
// 	rate_t_rp *modCurr;
// } OccupancyVectorIters;




struct ProducerElement {
	ElementInfo
	rate_t_rp decayConstant;
	struct {
		rate_t_rp productionConstant;
		pmbPtr produces;
		OccupancyVector occupancies;
		BoundElementArray boundelts;
		ProducerRxns reactions;
	};
	
};

struct ModulatorElement {
	ElementInfo
	rate_t_rp decayConstant;
	BoundElementPtrArray boundelts;
	effect_t_rp productionEffect;
	effect_t_rp decayEffect;
	Reaction * selfDecay;
	ReactionPtrArray bindings;
};

#ifndef DEBUG
inline pmbPtr pmbFromProducer(Producer * p){
	return (pmbPtr) {.producer = p, .species = p->species};
}

inline pmbPtr pmbFromModulator(Modulator * m){
	return (pmbPtr) {.modulator = m, .species = m->species};
}

inline ulong_type getElementID(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->id;
	else if(isModulator(toGet.species)) return toGet.modulator->id;
	else if(isBound(toGet.species)) return toGet.bound->id;
	else return ulong_type_max;
}

inline ulong_type getElementQty(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->qty;
	else if(isModulator(toGet.species)) return toGet.modulator->qty;
	else if(isBound(toGet.species)) return toGet.bound->qty;
	else return ulong_type_max;
}
#else
pmbPtr pmbFromProducer(Producer * p);
pmbPtr pmbFromModulator(Modulator * m);

ulong_type getElementID(pmbPtr toGet);
ulong_type getElementQty(pmbPtr toGet);
#endif

#endif