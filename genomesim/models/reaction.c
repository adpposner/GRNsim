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
		case MESS_TRANSCRIPTION: if(specLeft) assert(r->src->species & specLeft);
				if(specRight) assert(r->prod.species & specRight);break;
		case MICRO_TRANSCRIPTION: if(specLeft) assert(r->src->species  & specLeft);
		if(specRight) assert(r->prod.species & specRight);break;
		case TF_TRANSLATION: if(specLeft) assert(r->src->species  & specLeft);
		if(specRight) assert(r->prod.species & specRight);break;
		case MESS_DECAY: if(specLeft) assert(r->toDecay.species & specLeft);
							break;
		case MICRO_DECAY: if(specLeft) assert(r->toDecay.species  & specLeft);
							break;
		case TF_DECAY: if(specLeft) assert(r->toDecay.species  & specLeft);
							break;
		case MESSMIR_BINDING: if(specLeft) assert(r->left->species  & specLeft);
							if(specRight) assert(r->right->species & specRight);break;
		case TFDNA_BINDING: if(specLeft) assert(r->left->species  & specLeft);
		if(specRight) assert(r->right->species & specRight);break;
		case MESSMIR_UNBINDING: if(specLeft) assert(r->left->species  & specLeft);
		if(specRight) assert(r->right->species & specRight);break;
		case TFDNA_UNBINDING: if(specLeft) assert(r->left->species  & specLeft);
		if(specRight) assert(r->right->species & specRight);break;
	
		default: assert(0);
	}
}

void assertMessTxC(Reaction*r){assertRxnProps(r,MESS_TRANSCRIPTION,CODING,MESSENGER);}
void assertMirTxC(Reaction*r){assertRxnProps(r,MICRO_TRANSCRIPTION,NONCODING,MICRO);}
void assertAnyTxC(Reaction *r){assertRxnProps(r, ANY_TRANSCRIPTION, DNA, (MESSENGER | MICRO));}
void assertTxL(Reaction *r){assertRxnProps(r, TF_TRANSLATION, MESSENGER, PROTEIN);}
void assertMessDecay(Reaction *r){ assertRxnProps(r, MESS_DECAY, MESSENGER, NOSPECIES);}
void assertMirDecay(Reaction *r){ assertRxnProps(r, MICRO_DECAY, MICRO, NOSPECIES);}
void assertTFDecay(Reaction *r){ assertRxnProps(r, TF_DECAY, PROTEIN, NOSPECIES);}
void assertModulatorDecay(Reaction *r){assertRxnProps(r, MOD_DECAY, MICRO | PROTEIN, NOSPECIES);}
void assertBindingTFDNA(Reaction *r){assertRxnProps(r, TFDNA_BINDING, DNA, PROTEIN);}
void assertBindingMessMir(Reaction *r){assertRxnProps(r, MESSMIR_BINDING, MESSENGER, MICRO);}
void assertUnbindingTFDNA(Reaction *r){assertRxnProps(r, TFDNA_UNBINDING, DNA, PROTEIN);}
void assertUnbindingMessMir(Reaction *r){assertRxnProps(r, MESSMIR_UNBINDING, MESSENGER, MICRO);}


//  dependencyPack getDependenciesForProduction(Reaction *r){
// 	//if it is a mess production, then we fill in self
// 	//leftDeps becomes the rxns for the messenger
// 	//rightdeps is gone

// 	if(r->prod.species & MESSENGER)
// 		return (dependencyPack) {.self=r,.leftDeps = r->prod.producer->reactions};
// 	//If it is a modulator production, then we fill in self
// 	//leftDeps is now gone, and rightDeps is in play
// 	else
// 		return (dependencyPack) {.self = r, .leftDeps = {0}, .rightDeps = getReactionPtrArrayIters(&r->prod.modulator->bindings)};
// }

//  dependencyPack getDependenciesForDecay(decayReaction d){
// 	//here it is a little tricky - if an mRNA decays, we want to update its production,
// 	//decay,bindings, but do the initial set of unbindings
// 	if(d.toDecay.species & MESSENGER)
// 		return (dependencyPack) {.self=NULL,.leftDeps = d.toDecay.producer->reactions,.rightDeps={0}};
// 	else
// 		return (dependencyPack) {.self=NULL,.leftDeps={0},.rightDeps = getReactionPtrArrayIters(&d.toDecay.modulator->bindings)};
// 	//if a modulator decays, then we want to just update its rxns
// }

//  dependencyPack getDependenciesForBinding(bindingReaction b){
// 	//Here we know what's on the left and right, so this is fairly straightforward
// 	return (dependencyPack) {.self = NULL,.leftDeps = b.left->}
// }
//  dependencyPack getDependenciesForUnbinding(bindingReaction ub);


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
		case MESS_TRANSCRIPTION: return "MESS_TRANSCRIPTION";
case MICRO_TRANSCRIPTION: return "MICRO_TRANSCRIPTION";
case TF_TRANSLATION: return "TF_TRANSLATION";
case MESS_DECAY: return "MESS_DECAY";
case MICRO_DECAY: return "MICRO_DECAY";
case TF_DECAY: return "TF_DECAY";
case MESSMIR_BINDING: return "MESSMIR_BINDING";
case TFCODING_BINDING: return "TFCODING_BINDING";
case TFNONCODING_BINDING: return "TFNONCODING_BINDING";
case MESSMIR_UNBINDING: return "MESSMIR_UNBINDING";
case TFCODING_UNBINDING: return "TFCODING_UNBINDING";
case TFNONCODING_UNBINDING: return "TFNONCODING_UNBINDING";
	default: return "INVALID REACTION TYPE";
	}
}

static ReactionArray reactions;
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
	toRet.baseRate = src->productionConstant;
	toRet.src = src;
	if (src->species == CODING)
		toRet.rxn_type = MESS_TRANSCRIPTION;
	else if (src->species == NONCODING)
		toRet.rxn_type = MICRO_TRANSCRIPTION;
	else if (src->species == MESSENGER)
		toRet.rxn_type = TF_TRANSLATION;
	else
	{
		assert(0);
	}
	toRet.prod = src->produces;
	//toRet.dependencies = emptyPtrArray();
	//Don't init dependencies yet

	return toRet;
}

Reaction initBaseDecayReaction(pmbPtr * src) {
	Reaction toRet;
	assert(src->species & (MESSENGER | MICRO | PROTEIN));
	toRet.id = currRxnIndex++;
	switch (src->species)
	{
		case MESSENGER:
		toRet.rxn_type = MESS_DECAY;
			break;
		case MICRO:
		toRet.rxn_type = MICRO_DECAY;
			break;
		case PROTEIN:
		toRet.rxn_type = TF_DECAY;
			break;
	
		default: assert(0);
			break;
	}
	if (src->species & MESSENGER)
		toRet.baseRate = src->producer->decayConstant;
	else toRet.baseRate = src->modulator->decayConstant;
	toRet.toDecay = *src;
	//toRet.dependencies = emptyPtrArray();
	return toRet;
}


//Base association/dissociation rates are based on symmetry about 1 for the time being
 rate_t_rp assocConst(rate_t_rp aff) {
	rate_t_rp k; k = 2.0 / (1.0 + aff); return 2.0 - k;
}
 rate_t_rp dissocConst(rate_t_rp aff) {
	rate_t_rp k; k = 2.0 / (1.0 + aff); return k;
}



Reaction initBindingReaction(BoundElement * boundElt) {
	Reaction toRet;
	toRet.id  = currRxnIndex++;
	switch(boundElt->left->species){
		case MESSENGER:
			toRet.rxn_type = MESSMIR_BINDING;
		break;
		case CODING:
			toRet.rxn_type = TFCODING_BINDING;
		break;
		case NONCODING:
			toRet.rxn_type = TFNONCODING_BINDING;
		break;
		default:
		assert(0);
	}

	toRet.baseRate = boundElt->assocConstant;
	toRet.left = boundElt->left;
	toRet.right = boundElt->right;
	toRet.target = boundElt;
	//toRet.dependencies = emptyPtrArray();
	return toRet;
}

Reaction initUnbindingReaction(BoundElement * boundElt) {
	Reaction toRet;
	toRet.id  = currRxnIndex++;
	switch(boundElt->left->species){
		case MESSENGER:
			toRet.rxn_type = MESSMIR_UNBINDING;
		break;
		case CODING:
			toRet.rxn_type = TFCODING_UNBINDING;
		break;
		case NONCODING:
			toRet.rxn_type = TFNONCODING_UNBINDING;
		break;
		default:
		assert(0);
	}
	toRet.baseRate = boundElt->dissocConstant;
	toRet.left = boundElt->left;
	toRet.right = boundElt->right;
	toRet.target = boundElt;
	//toRet.dependencies = emptyPtrArray();
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
		case ANY_PRODUCTION: return (isProducer(pmb.species)) ? 1 : 0;
		case ANY_DECAY: return (canDecay(pmb.species)) ? 1 : 0;
		case ANY_UNBINDING:
		case ANY_BINDING: return isBound(pmb.species) ? 1 : 0;
		default: perror("InvalidNumReactionTypeForReactant"); exit(-2313);
		}
	} else return 0;
}

static  unsigned totalReactionsForReactant(const pmbPtr pmb) {
	reaction_t rxntypes[4] = {ANY_PRODUCTION, ANY_DECAY, ANY_BINDING, ANY_UNBINDING};
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
	pmbPtr currElem = emptyPmbPtr();
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
//and unbindings which will be init'd IN THAT ORDER
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
	if (numReactionTypeForReactant(pmb, ANY_PRODUCTION)) *it.curr++ = initProductionReaction(pmb.producer);
	if (numReactionTypeForReactant(pmb, ANY_DECAY)) *it.curr++ = initBaseDecayReaction((pmbPtr *)&pmb);
	it.curr = initReactionsForBoundElementArray(&pmb.producer->boundelts, it.curr);

	assert((it.curr - currPos) == totalRxnsToInit);
	assert(it.curr <= it.end);
	return it.curr;
}

Reaction * initReactionsFromModulator(const pmbPtr pmb, ReactionArray * arr, Reaction * currPos) {
	assert(isModulator(pmb.species));
	ulong_type totalRxnsToInit = numReactionTypeForReactant(pmb, ANY_DECAY);
	if(!totalRxnsToInit) return currPos;
	assert(totalRxnsToInit);
	ReactionArrayIters it = getReactionArrayIters(arr);
	assert((currPos - it.start) >= 0);
	assert((it.end - currPos) > 0);
	it.curr = currPos;
	if (numReactionTypeForReactant(pmb, ANY_DECAY)) *it.curr++ = initBaseDecayReaction((pmbPtr *)&pmb);
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
	if ((ra & ANY_PRODUCTION) || (ra & ANY_UNBINDING)) dca = 1;
	if ((rb & ANY_PRODUCTION) || (rb & ANY_UNBINDING)) dcb = 1;
	if ((ra & ANY_DECAY) && ((*pa)->toDecay.species & MESSENGER)) dca = 1;
	if ((rb & ANY_DECAY) && ((*pb)->toDecay.species & MESSENGER)) dcb = 1;
	if (dca || dcb) return dca - dcb;
	ulong_type ia = (ra & ANY_DECAY) ? (*pa)->toDecay.modulator->id :
	                (*pa)->right->id;
	ulong_type ib = (rb & ANY_DECAY) ? (*pb)->toDecay.modulator->id :
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
        case MESS_TRANSCRIPTION:
        case MICRO_TRANSCRIPTION:
        case TF_TRANSLATION: ia = pa->src->id; rankA = 0; break;
        case MESS_DECAY: ia=pa->toDecay.producer->id; rankA=2;break;
        case TF_DECAY:
        case MICRO_DECAY: ia = pa->toDecay.modulator->id; dca = 1; rankA=2;break;
        case MESSMIR_BINDING:
        case TFCODING_BINDING:
        case TFNONCODING_BINDING: ia = pa->left->id; rankA = 1; break;
        case MESSMIR_UNBINDING:
        case TFCODING_UNBINDING:
        case TFNONCODING_UNBINDING: ia = pa->left->id; rankA = 3; break;
	    default: exit(-134); break;
	}
	switch (rb) {
	 case MESS_TRANSCRIPTION:
        case MICRO_TRANSCRIPTION:
        case TF_TRANSLATION: ib = pb->src->id; rankB = 0; break;
	    case MESS_DECAY: ib = pb->toDecay.producer->id; rankB = 2;break;
		case TF_DECAY:
        case MICRO_DECAY:  ib = pb->toDecay.producer->id; dcb = 1;
		    rankB = 2; break;
	    case MESSMIR_BINDING:
        case TFCODING_BINDING:
        case TFNONCODING_BINDING: ib = pb->left->id; rankB = 1; break;
	    case MESSMIR_UNBINDING:
        case TFCODING_UNBINDING:
        case TFNONCODING_UNBINDING: ib = pb->left->id; rankB = 3; break;
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
	assert(binding->rxn_type & ANY_BINDING);
	if (unbindingStart == unbindingEnd)
		return;
	const Reaction * unbindingCurr;
	for (unbindingCurr = unbindingStart; unbindingCurr != unbindingEnd; unbindingCurr++) {
		assert(unbindingCurr->rxn_type & ANY_UNBINDING);
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
		assert((*rIt.curr)->rxn_type & MOD_DECAY);
		mIt.curr->selfDecay = *rIt.curr++;
		startPos = rIt.curr;
		ReactionPtr * startPos = rIt.curr;
		while ((*rIt.curr)->rxn_type & ANY_BINDING) {
			nBindingReactions++;
			rIt.curr++;
		}
		endPos = rIt.curr;
		mIt.curr->bindings = ReactionPtrArray_alloc(nBindingReactions);
		ReactionPtrArrayIters rpIt = getReactionPtrArrayIters(&mIt.curr->bindings);
		rpIt.curr = rpIt.start;
		for (rIt.curr = startPos; rIt.curr != endPos; rIt.curr++) {
			*rpIt.curr++ = *rIt.curr;
			assert((*rIt.curr)->rxn_type & ANY_BINDING);
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
		assert(rIt.curr->rxn_type & ANY_PRODUCTION);
		pIt.curr->reactions.production = rIt.curr++;
		if (rIt.curr->rxn_type & ANY_BINDING) {
			pIt.curr->reactions.binding = rIt.curr++;
			assert(rIt.curr->left == pIt.curr);
		}
		while (rIt.curr->rxn_type & ANY_BINDING) {
			assert(rIt.curr->left == pIt.curr);
			rIt.curr++;
		}
		if (pIt.curr->species & MESSENGER) {
			assert(rIt.curr->left == pIt.curr);
			assert(rIt.curr->rxn_type & ANY_DECAY);
			pIt.curr->reactions.decay = rIt.curr++;

		}
		if (rIt.curr->rxn_type & ANY_UNBINDING) {
			assert(rIt.curr->left == pIt.curr);
			pIt.curr->reactions.unbinding = rIt.curr++;
		}
		while (rIt.curr->rxn_type & ANY_UNBINDING) {
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
///BEGIN REACTION CONNECTION ROUTINES ///////

// static  unsigned char isReactantIn(Reaction * r, BasicGeneticElement *bge){

// 	switch(r->rxn_type){
// 		case PRODUCTION: return (r->src->id == bge->id) ? 1 : 0;
// 		case DECAY: return (r->toDecay->id == bge->id) ?  1 : 0;
// 		case BINDING:	if (r->left->id == bge->id) return 1;
// 						else if (r->right->id == bge->id) return 1;
// 						return 0;
// 		case UNBINDING: return (r->target->id == bge->id) ? 1 : 0;
// 	}
// }

// static unsigned char * isSet=NULL;
// void markReactionsFor(BasicGeneticElement * bge,unsigned rxntypes) {

// 	if(!isSet) isSet = calloc(reactions.length,sizeof(*isSet));
// 	ReactionArrayIters it = getReactionArrayIters(&reactions);
// 	unsigned ind;
// 	for(it.curr=it.start,ind=0;it.curr!=it.end;it.curr++,ind++)
// 		if(!(rxntypes & it.curr->rxn_type)) continue;
// 		else isSet[ind]=isSet[ind]+isReactantIn(it.curr, bge);
// }

// ReactionPtrArray exportRxnPtrArray() {
// 	unsigned ct,i;

// 	for(ct=0, i = 0; i<reactions.length;i++)
// 		if(isSet[i]) ct++;
// 	ReactionPtrArray toRet;

// 	toRet = ReactionPtrArray_alloc(ct);
// 	ReactionPtrArrayIters itTarget = getReactionPtrArrayIters(&toRet);
// 	itTarget.curr=itTarget.start;
// 	for(i=0;i<reactions.length;i++)
// 		if (isSet[i]) *itTarget.curr++ = &reactions.data[i];
// 	assert(itTarget.curr==itTarget.end);
// 	memset(isSet, 0, reactions.length);
// 	return toRet;
// }

// ReactionPtrArray getProductionDependencies(Reaction * r) {

// 	markReactionsFor(r->src, PRODUCTION); //self
// 	markReactionsFor(r->prod, PRODUCTION | DECAY | BINDING);
// 	return exportRxnPtrArray();
// }

// ReactionPtrArray getDecayDependencies(Reaction * r) {
// 	BasicGeneticElement * tD;
// 	tD = r->toDecay;
// 	markReactionsFor(tD, DECAY | PRODUCTION | BINDING);
// 	//attach unbindings to decay
// 	if(tD->species & (DNA | MESSENGER)){
// 		BasicGeneticElementArrayIters it = getBasicGeneticElementArrayIters(&tD->base.base.boundelts.elements);
// 		for(it.curr=it.start;it.curr!=it.end;it.curr++)
// 			markReactionsFor(it.curr, UNBINDING);
// 	}

// 	return exportRxnPtrArray();
// }

// ReactionPtrArray getBindingUnbindingDependencies(Reaction *r) {
// 	markReactionsFor(r->left, PRODUCTION);
// 	markReactionsFor(r->right, DECAY | BINDING);
// 	markReactionsFor(r->target, UNBINDING);
// 	return exportRxnPtrArray();
// }

// void setDependencies(Reaction *r) {
// 	switch (r->rxn_type){
// 		case PRODUCTION: r->dependencies = getProductionDependencies(r);break;
// 		case DECAY: r->dependencies = getDecayDependencies(r);break;
// 		case BINDING:
// 		case UNBINDING: r->dependencies = getBindingUnbindingDependencies(r);break;
// 	}
// }

// void generateReactionDependencies(ReactionArray *rxns) {
// 	ReactionArrayIters it = getReactionArrayIters(rxns);
// 	for(it.curr=it.start;it.curr!=it.end;it.curr++){
// 		setDependencies(it.curr);
// 	}
// }

/////END DEPENDENCY GENERATION /////////////////////
///BEGIN UTILITIES	///////////////////////////
void printRxn(FILE * fh, const Reaction * rxn, const int wDeps) {
#define PRINTRXNBUFSIZE	1000
	char buf[PRINTRXNBUFSIZE] = {0};
	memset(buf, 0, PRINTRXNBUFSIZE);
	int cx = 0;
	int n = 0;

	productionReaction p;
	decayReaction d;
	bindingReaction b;

	n += snprintf(buf + cx, PRINTRXNBUFSIZE - cx, "ID = %u\t", rxn->id); if (n >= 0)cx += n;
    reaction_t rt = rxn->rxn_type;
	if(rt & ANY_PRODUCTION){cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s: ", event_name(rxn->rxn_type));

		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "SRCTYPE = %s SRC ID = %zu, PROD ID = %zu SRCQTY= %zu, PRODQTY = %zu", species_names(rxn->src->species),
		               rxn->src->id, getElementID(rxn->prod), rxn->src->qty, getElementQty(rxn->prod));}
	else if(rt & ANY_DECAY)
	{cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s: ", event_name(rxn->rxn_type));
		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "SPTYPE = %s ID = %zu QTY = %zu ",
		               species_names(rxn->toDecay.species), getElementID(rxn->toDecay), getElementQty(rxn->toDecay));
    }else if (rt & ( ANY_BINDING | ANY_UNBINDING))
	{ cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "%s", event_name(rxn->rxn_type));
		cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, "LEFT: %zu, RIGHT:%s, %zu, TARG: %zu QTY = %zu ", rxn->left->id, (rxn->right->species == PROTEIN) ? "PROT" : "MIR", rxn->right->id, rxn->target->id, rxn->target->qty);
    }else fprintf(fh, "UNKNOWN REACTION TYPE %s:%d TYPE = %d", __FILE__, __LINE__, rxn->rxn_type); //exit(-5);
	
	cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, " TTE: %f RATE: %f", rxn->timeToEvent, rxn->reactionRate);
	cx += snprintf(buf + cx, PRINTRXNBUFSIZE - 1 - cx, " DEPS: ");
	//ReactionPtrArrayIters it = getReactionPtrArrayIters(&rxn->dependencies);
	// if (!wDeps)
	// {	for(it.curr=it.start;it.curr!=it.end;it.curr++)
	// 		{n+=snprintf(buf+cx,PRINTRXNBUFSIZE-1-cx,"%u, ",(**it.curr).id);if(n>=0)cx+=n;}
	// 		fprintf(fh,"%s\n",buf);
	//    if (wDeps) {fprintf(fh,"%s\n",buf);
	// for(it.curr=it.start;it.curr!=it.end;it.curr++)
	// 	printRxn(fh,*it.curr,0);}
	//else
	fprintf(fh, "\t%s\n", buf);


}

void printAllRxns(FILE * fh, ReactionArray *r) {
	ReactionArrayIters it = getReactionArrayIters(r);
	ARRAY_TYPE_FOREACH(it) printRxn(fh, it.curr, 0);
}

