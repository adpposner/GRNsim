#include "../include/globals.h"
#include "../include/models/basicgeneticelement.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/parameters.h"
#include "../include/randomnumbers.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../include/globals.h"
#include "../include/models/reaction.h"
#include <math.h>
#include <limits.h>




char isProducer(species_t species) {
	return (species & (CODING | NONCODING | MESSENGER)) ? 1 : 0;
}

char isModulator(species_t species) {
	return (species & (MICRO | PROTEIN)) ? 1 : 0;
}

char isBound(species_t species) {
	return (species & (MESSMIR | TFDNA)) ? 1 : 0;
}


static const occ_bits_rp_t nullElem = {.x = {1<<0,0,0,0}};

char canDecay(species_t species) {
	#ifdef BOUNDDECAY
	#define BDCANDECAY	MESSMIR | TFDNA
	#else
	#define BDCANDECAY	0
	#endif
	return (species & (MESSENGER | MICRO | PROTEIN | BDCANDECAY)) ? 1 : 0;
}

ProducerRxns nullProducerRxns(){
	return (ProducerRxns){.production = NULL,.decay = NULL,.binding = NULL,
		.unbinding=NULL,.end=NULL};
}


pmbPtr pmbFromBound(BoundElement * b){
	return (pmbPtr) {.bound = b, .species = b->species};
}


#ifdef DEBUG
pmbPtr pmbFromProducer(Producer * p){
	return (pmbPtr) {.producer = p, .species = p->species};
}

pmbPtr pmbFromModulator(Modulator * m){
	return (pmbPtr) {.modulator = m, .species = m->species};
}

ulong_type getElementID(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->id;
	else if(isModulator(toGet.species)) return toGet.modulator->id;
	else if(isBound(toGet.species)) return toGet.bound->id;
	else return ulong_type_max;
}

ulong_type getElementQty(pmbPtr toGet){
	if (isProducer(toGet.species)) return toGet.producer->qty;
	else if(isModulator(toGet.species)) return toGet.modulator->qty;
	else if(isBound(toGet.species)) return toGet.bound->qty;
	else return ulong_type_max;
}
#endif
OccupancyVector emptyOccupancyVector() {
	return (OccupancyVector) {.len = 0, .data = NULL,.decaySum=0.0,.prodSum=0.0};
}

void OccupancyVector_free(OccupancyVector * toFree) {
	if (toFree->data)
		free(toFree->data);
	toFree->len = 0;
	toFree->data = NULL;
}

static OccupancyVector OccupancyVector_alloc(unsigned maxLen) {
	OccupancyVector toRet = emptyOccupancyVector();
	toRet.len = maxLen;
	if (maxLen) {
		toRet.data = malloc(sizeof(*toRet.data) * maxLen);
		memset(toRet.data, 0, toRet.len * sizeof(*toRet.data));
	}
	return toRet;
}

static void initOccupancyVector(Producer *p, ulong_type initialQty, unsigned maxVecLen) {
	assert(initialQty <= maxVecLen);
	p->occupancies = OccupancyVector_alloc(maxVecLen);
	unsigned i;
	p->occupancies.prodSum=0.0;
	p->occupancies.decaySum=0.0;
	for (i = 0; i < initialQty; i++) {
		p->occupancies.data[i].bits = nullElem;
		p->occupancies.data[i].prodMod = 1.0;
		p->occupancies.data[i].decayMod = 1.0;
		p->occupancies.prodSum++;
		p->occupancies.decaySum++;
	}
}




static void BoundElement_free(BoundElement *toFree) {}
static void BoundElementPtr_free(BoundElementPtr *toFree) {}
static void Producer_free(Producer * toFree) {
	OccupancyVector_free(&toFree->occupancies);
}
static void Modulator_free(Modulator *toFree) {

	BoundElementPtrArray_free(&toFree->boundelts);
	ReactionPtrArray_free(&toFree->bindings);
}
static void ProducerPtr_free(ProducerPtr * p) {}
static void ModulatorPtr_free(ModulatorPtr * m) {}
static void pmbPtr_free(pmbPtr * tpFree) {}

DEFINEBASICARRAYTYPEPRIMITIVE(BoundElement, BoundElement)
DEFINEBASICARRAYTYPEPRIMITIVE(BoundElementPtr, BoundElementPtr)
DEFINEBASICARRAYTYPEPRIMITIVE(Producer, Producer)
DEFINEBASICARRAYTYPEPRIMITIVE(Modulator, Modulator)
DEFINEBASICARRAYTYPEPRIMITIVE(ProducerPtr, ProducerPtr)
DEFINEBASICARRAYTYPEPRIMITIVE(ModulatorPtr, ModulatorPtr)
DEFINEBASICARRAYTYPEPRIMITIVE(pmbPtr, pmbPtr)

static ulong_type nextID() {
	static ulong_type ctr = 0;
	return ctr++;
}

static ulong_type nextName() {
	static ulong_type tfn = 0;
	return tfn++;
}



pmbPtr emptyPmbPtr() {return (pmbPtr) {.producer = NULL, .species = UNDEFINED};}

static  pmbPtr ptrForMess(Producer * messRNA) {return (pmbPtr) {.producer = messRNA, .species = MESSENGER};}
static  pmbPtr ptrForMicro(Modulator * micro) {return (pmbPtr) {.modulator = micro, .species = MICRO};}
static  pmbPtr ptrForProt(Modulator * prot) {return (pmbPtr) {.modulator = prot, .species = PROTEIN};}
static  pmbPtr ptrForPmb(void * targ, species_t spec) {
	switch (spec) {
	case MESSENGER: return (pmbPtr) {.producer = (Producer *) targ, .species = spec};
	case PROTEIN:
	case MICRO: return (pmbPtr) {.modulator = (Modulator *) targ, .species = spec};
	case MESSMIR:
	case TFDNA: return (pmbPtr) {.bound = (BoundElement *) targ, .species = spec};
	default:
		exit(-1234125);
	}
}

static  ulong_type getIDFromPMBPtr(pmbPtr * p) {
	if (p->species & (CODING | NONCODING | MESSENGER))
		return p->producer->id;
	else if (p->species & (MICRO | PROTEIN))
		return p->modulator->id;
	else if (p->species & (TFDNA | MESSMIR))
		return p->bound->id;
	else assert(0);
	return 0;
}

ulong_type getProducesID(Producer * p) {return getIDFromPMBPtr(&p->produces);}






void setDefaultProducerQuantity(Producer * p) {
	setProducerQty(p, getInitialQty(p->species));
}
void setDefaultModulatorQuantity(Modulator *m) {
	m->qty = getInitialQty(m->species);
}

void setDefaultBoundQuantity(BoundElement * b) {
	b->qty = getInitialQty(b->species);
}

void setProducerQty(Producer * p, unsigned initialQty) {
	p->qty = initialQty;
	unsigned maxsize;
	OccupancyVector_free(&p->occupancies);
	switch (p->species) {
	case CODING: p->qty = maxsize = CODING_QTY;
		if (initialQty != CODING_QTY)
			fprintf(stderr, "Requested %u but defaulting to %u for coding", initialQty, CODING_QTY);
		break;
	case NONCODING: p->qty = maxsize = NONCODING_QTY;
		if (initialQty != NONCODING_QTY)
			fprintf(stderr, "Requested %u but defaulting to %u for noncoding", initialQty, NONCODING_QTY); break;
	case MESSENGER: p->qty = initialQty;
		if (initialQty > MAX_MRNA_QTY) {fprintf(stderr, "Too many miRNAs - IC = %u, max = %u", initialQty, MAX_MRNA_QTY); exit(-690782);}
		maxsize = MAX_MRNA_QTY; break;
	default: exit(-9000);
	}
	initOccupancyVector(p, p->qty, maxsize);
}



static Producer initDNA(char * name, ulong_type id, species_t spec, void * product,
                        rate_t_rp productionConstant, BoundElementArray * conns, unsigned char enabled) {
	Producer toRet;
	toRet.id = id;
	toRet.species = spec;
	toRet.qty = 0;
	toRet.name[0] = '\0';
	if (name)
		snprintf(toRet.name, MAXGENENAMELENGTH, "%s", name);
	toRet.produces = ptrForPmb(product, (spec == CODING) ? MESSENGER : MICRO);
	toRet.productionConstant = productionConstant;
	toRet.decayConstant = 0;
	toRet.boundelts = (conns) ? *conns : emptyBoundElementArray();
	toRet.isEnabled = enabled;
	toRet.occupancies = emptyOccupancyVector();
	return toRet;
}

Producer initCoding(char * name, ulong_type id, Producer * messRna,
                    rate_t_rp * productionConstant, BoundElementArray * conns, unsigned char enabled) {
	return initDNA(name, id, CODING, messRna, messengerProductionRate(productionConstant), conns, enabled);
}

Producer initNoncoding(char * name, ulong_type id, Modulator * microRna,
                       rate_t_rp * productionConstant, BoundElementArray * conns, unsigned char enabled) {
	return initDNA(name, id, NONCODING, microRna, microProductionRate(productionConstant), conns, enabled);
}

Producer initMessenger(char * name, ulong_type id, Modulator * tfprotein,
                       rate_t_rp * productionConstant, rate_t_rp * decayConstant, BoundElementArray * conns, unsigned char enabled) {
	Producer toRet;
	toRet.id = id;
	toRet.species = MESSENGER;
	toRet.qty = 0;
	toRet.name[0] = '\0';
	if (name)
		snprintf(toRet.name, MAXGENENAMELENGTH, "%s", name);
	toRet.produces = ptrForProt(tfprotein);
	toRet.productionConstant = proteinProductionRate(productionConstant);
	toRet.decayConstant = messengerDecayRate(decayConstant);
	toRet.isEnabled = enabled;
	toRet.boundelts = (conns) ? *conns : emptyBoundElementArray();
	toRet.occupancies = emptyOccupancyVector();
	return toRet;
}

Modulator initModulator(char *name, ulong_type id, species_t spec, rate_t_rp decayConstant,
                        BoundElementPtrArray *conns, effect_t_rp prodEffect, effect_t_rp decayEffect, unsigned char enabled) {
	Modulator toRet;
	toRet.id = id;
	toRet.species = spec;
	toRet.qty = 0;
	toRet.name[0] = '\0';
	if (name)
		snprintf(toRet.name, MAXGENENAMELENGTH, "%s", name);
	toRet.productionEffect = prodEffect;
	toRet.decayEffect = decayEffect;
	toRet.decayConstant = decayConstant;
	toRet.boundelts = (conns) ? *conns : emptyBoundElementPtrArray();
	toRet.bindings = emptyReactionPtrArray();
	toRet.isEnabled = enabled;
	return toRet;
}


Modulator initMiRNA(char * name, ulong_type id,
                    rate_t_rp * decayConstant, effect_t_rp * prodEffect,effect_t_rp * decayEffect,
					 BoundElementPtrArray * conns, unsigned char enabled) {
	return initModulator(name, id, MICRO, microDecayRate(decayConstant), conns,miRProdEffectStrength(prodEffect),miRDecayEffectStrength(decayEffect),enabled);
}

Modulator initTF(char * name, ulong_type id, rate_t_rp * decayConstant, effect_t_rp * prodEffect,
                 BoundElementPtrArray * conns,unsigned char enabled) {
	return initModulator(name, id, PROTEIN, proteinDecayRate(decayConstant), conns,TFEffectStrength(prodEffect),0.0, enabled);
}

BoundElement initBound(char * name, ulong_type id, Producer * left, species_t spec,
                       Modulator * right, rate_t_rp assocConst, rate_t_rp dissocConst,const rate_t_rp *decayConst, 
					   effect_t_rp prodEffectStrength,
					   effect_t_rp decayEffectStrength) {
	BoundElement toRet;
	toRet.id = id;
	toRet.species = spec;
	toRet.qty = 0;
	if (name) snprintf(toRet.name, MAXGENENAMELENGTH, "%s", name);
	else snprintf(toRet.name, MAXGENENAMELENGTH, "%s-%s", left->name, right->name);
	toRet.left = left;
	toRet.right = right;
	toRet.assocConstant = assocConst;
	toRet.dissocConstant = dissocConst;
	toRet.prodEffectStrength = prodEffectStrength;
	toRet.decayEffectStrength = decayEffectStrength;
	toRet.producerArrayPos = 0;
	toRet.unbinding = NULL;
#ifdef BOUNDDECAY
	toRet.decayConstant = (decayConst) ? *decayConst : 0.0;
#endif
	return toRet;
}

BoundElement initMessMir(char * name, ulong_type id, Producer * left,
                         Modulator * right, rate_t_rp *assocConst, rate_t_rp *dissocConst, rate_t_rp *decayConst, 
						 const effect_t_rp *prodEffectStrength,const effect_t_rp * decayEffectStrength) {
	assert(left->species == MESSENGER); assert(right->species == MICRO);
	return initBound(name, id, left, MESSMIR, right, miRAssocConst(assocConst), miRDissocConst(dissocConst), decayConst, 
	miRProdEffectStrength(prodEffectStrength),miRDecayEffectStrength(decayEffectStrength));

}

BoundElement initTFDNA(char * name, ulong_type id, Producer * left,
                       Modulator * right,const rate_t_rp *assocConst,const rate_t_rp *dissocConst, 
					   const rate_t_rp *decayConst,const effect_t_rp *effectStrength) {
	assert(left->species & DNA); assert(right->species == PROTEIN);
	return initBound(name, id, left, TFDNA, right, TFAssocConst(assocConst), TFDissocConst(dissocConst), decayConst, 
	TFEffectStrength(effectStrength),0.0);
}

static BoundElement createBasicBound(Producer * l, Modulator * r)
{
	BoundElement toRet;

	species_t rightSpecies = r->species;
	species_t leftSpecies = l->species;
	if (rightSpecies & PROTEIN) {
		assert(leftSpecies & DNA);
		return initTFDNA(NULL, nextID(), l, r, NULL, NULL, NULL, &r->productionEffect);
	}
	else if (rightSpecies == MICRO) {
		assert(leftSpecies == MESSENGER);
		return initMessMir(NULL, nextID(), l, r, NULL, NULL, NULL, &r->productionEffect,&r->decayEffect);
	} else {
       fprintf(stderr,"%d\t%d",leftSpecies,rightSpecies);
        perror("Failure in createBasicBound"); exit(-1023);}

}

BoundElement * initConnection(pmbPtr *src, pmbPtr * targ, BoundElement * pos) {
	if (src->species & PRODUCER_SPECIES) {
		assert(targ->species & MODULATOR_SPECIES);
		*pos++ = createBasicBound(src->producer, targ->modulator);
	} else if (src->species & MODULATOR_SPECIES) {
		assert(targ->species & PRODUCER_SPECIES);
		*pos++ = createBasicBound(targ->producer, src->modulator);
	} else
		exit(-1253);
	return pos;
}

void initBasicTF(char *name, Producer * dna, Producer * mess, Modulator * tf, unsigned char enabled) {
	char geneName[MAXGENENAMELENGTH];
	char messName[MAXGENENAMELENGTH];
	char tfName[MAXGENENAMELENGTH];
	if (name) {
		snprintf(geneName, MAXGENENAMELENGTH, "gene %s", name);
		snprintf(messName, MAXGENENAMELENGTH, "mess %s", name);
		snprintf(tfName, MAXGENENAMELENGTH, "prot %s", name);
	} else {
		ulong_type nName = nextName();
		snprintf(geneName, MAXGENENAMELENGTH, "gene %zu", nName);
		snprintf(messName, MAXGENENAMELENGTH, "mess %zu", nName);
		snprintf(tfName, MAXGENENAMELENGTH, "prot %zu", nName);
	}
	*dna = initCoding(geneName, nextID(), mess, NULL, NULL, enabled);
	*mess = initMessenger(messName, nextID(), tf, NULL, NULL, NULL, enabled);
	*tf = initTF(tfName,nextID(),NULL,NULL,NULL,enabled);
	
}

void initBasicMicro(char *name, Producer *dna, Modulator * miR, unsigned char enabled) {
	char ncgeneName[MAXGENENAMELENGTH];
	char miRName[MAXGENENAMELENGTH];
	if (name) {
		snprintf(ncgeneName, MAXGENENAMELENGTH, "ncgene %s", name);
		snprintf(miRName, MAXGENENAMELENGTH, "mir %s", name);

	} else {
		ulong_type nName = nextName();
		snprintf(ncgeneName, MAXGENENAMELENGTH, "ncgene %zu", nName);
		snprintf(miRName, MAXGENENAMELENGTH, "mir %zu", nName);

	}
	*dna = initNoncoding(ncgeneName, nextID(), miR, NULL, NULL, enabled);
	*miR = initMiRNA(miRName,nextID(),NULL,NULL,NULL,NULL,enabled);
}

void Element_print(FILE * fp, pmbPtr elem) {
	char buf[1000];
	int n = 0;
	int cx = 0;

	Producer *p = NULL;
	Modulator *m = NULL;
	BoundElement *b = NULL;
	switch (elem.species) {
	case UNDEFINED: assert(0); break;
	case CODING:
	case NONCODING:
	case DNA:
	case MESSENGER:
		p = elem.producer;
		n += snprintf(buf, 1000, "%s\t%zu\t", species_names(elem.species), m->id); if (n > 0) cx += n;
		n += snprintf(buf + cx, 1000 - cx, "produces %zu\tqty=%zu\n", getProducesID(p), p->qty);
		break;
	case MICRO:
	case PROTEIN:
		m = elem.modulator;
		n += snprintf(buf, 1000, "%s\t%zu\t", species_names(elem.species), m->id); if (n > 0) cx += n;
		n += snprintf(buf + cx, 1000 - cx, "qty = %zu\n", m->qty);
		break;
	case MESSMIR:
	case TFDNA:
		b = elem.bound;
		n += snprintf(buf, 1000, "%s\t%zu\t", species_names(elem.species), b->id); if (n > 0) cx += n;
		n += snprintf(buf + cx, 1000 - cx, "left=%zu\tright=%zu\tqty=%zu\n", b->left->id,
		              b->right->id, b->qty); break;
	default: exit(-34782);
	}
	fprintf(fp, "%s", buf);

}

static int BoundElementCmpFunc(const void *a, const void *b) {
	BoundElement *pa = (BoundElement *)a;
	BoundElement *pb = (BoundElement *)b;
	ulong_type aLeft, aRight, bLeft, bRight;
	aLeft = pa->left->id;
	aRight = pa->right->id;
	bLeft = pb->left->id;
	bRight = pb->right->id;
	if (aLeft < bLeft) return -1;
	else if (bLeft < aLeft) return 1;
	else {
		if (aRight < bRight) return -1;
		else if (aRight > bRight) return 1;
		else assert(0);
	}
	return 0;
}

void sortBoundElements(BoundElementArray bd) {
	qsort(bd.data, bd.length, sizeof(*bd.data), BoundElementCmpFunc);
}


//7.20 now assigns arrayIndex to bound element, each should be unique now
void assignBoundElementsToProducer(Producer *p, const BoundElementArray *bd) {
	//assume sorted
	assert(p->boundelts.data == NULL);
	assert(p->boundelts.length == 0);
	BoundElementArrayIters it = getBoundElementArrayIters(bd);
	for (it.curr = it.start; it.curr != it.end; it.curr++) {
		if ((p->boundelts.length == 0) && (p == it.curr->left)) p->boundelts.data = it.curr;
		if (p == it.curr->left) {
			p->boundelts.length++;
			assert(it.curr->producerArrayPos == 0);
			it.curr->producerArrayPos = (unsigned char) p->boundelts.length;
			assert((p->boundelts.length < UCHAR_MAX) && (it.curr->producerArrayPos != 0));
		}
	}
	//perform check to ensure correctness - must have all be consecutive
#ifdef TESTASSERTIONS
	ulong_type lastID = 0;
	BoundElementArrayIters bIt = getBoundElementArrayIters(&p->boundelts);
	for (bIt.curr = bIt.start; bIt.curr != bIt.end; bIt.curr++) {
		assert(bIt.curr->left == p);
		assert(lastID < bIt.curr->right->id);
		lastID = bIt.curr->right->id;
	}
#endif
}

ulong_type countElementsWithRight(const Modulator * m, const BoundElementArray *bd) {
	ulong_type toRet = 0;
	BoundElementArrayIters it = getBoundElementArrayIters(bd);
	for (it.curr = it.start; it.curr != it.end; it.curr++)
		if (m == it.curr->right) toRet++;
	return toRet;
}

void initBoundElementPtrsForModulator(Modulator *m, const BoundElementArray *bd) {
	assert(m->boundelts.data == NULL);
	assert(m->boundelts.length == 0);
	m->boundelts.length = countElementsWithRight(m, bd);
	m->boundelts = BoundElementPtrArray_alloc(m->boundelts.length);
	BoundElementPtrArrayIters mIt = getBoundElementPtrArrayIters(&m->boundelts);
	BoundElementArrayIters it = getBoundElementArrayIters(bd);
	for (it.curr = it.start, mIt.curr = mIt.start; it.curr != it.end; it.curr++)
	{if (m == it.curr->right) {assert(m->id == it.curr->right->id); *mIt.curr++ = it.curr;}}
	assert(mIt.curr == mIt.end);

#ifdef TESTASSERTIONS

	ulong_type lastID = 0;
	BoundElementPtrArrayIters qIt = getBoundElementPtrArrayIters(&m->boundelts);
	if (qIt.curr != qIt.end) {
		lastID = (*qIt.curr)->left->id;
		qIt.curr++;
		for (; qIt.curr != qIt.end; qIt.curr++) {
			if (lastID >= (*qIt.curr)->left->id) {
				printf("%zu >= %zu, right = %zu, targ = %zu\n", lastID, (*qIt.curr)->left->id,
				       (*qIt.curr)->right->id, (*qIt.curr)->id);

				exit(-272);
			}
			lastID = (*qIt.curr)->left->id;
		}
	}
#endif
}



//Get degrees
ulong_type getIndegree_dna(const Producer *p) {
	assert(p->species & DNA);
	return p->boundelts.length;
}

ulong_type getIndegree_mess(const Producer *mess) {
	assert(mess->species == MESSENGER);
	BoundElementArrayIters bIt = getBoundElementArrayIters(&mess->boundelts);
	ulong_type toRet = 0;
	ARRAY_TYPE_FOREACH(bIt) if (bIt.curr->right->isEnabled) toRet++;
	return toRet;
}
//end get degrees


void calcMaxEffect(Producer * p,numeric_t_rp * maxProdEffect, numeric_t_rp * baseProdEffect, numeric_t_rp * maxDecayEffect,numeric_t_rp * baseDecayEffect){
	BoundElementArrayIters bIt = getBoundElementArrayIters(&p->boundelts);
	numeric_t_rp pE,dE,bP,bD;
	pE = dE = bP = bD = 1.0;
	ARRAY_TYPE_FOREACH(bIt){
		pE *= bIt.curr->prodEffectStrength;
		dE *= bIt.curr->decayEffectStrength;
	}
	if(p->species & DNA){
		pE = pow(pE,TXC_PROD_EXPONENT);
		pE = p->productionConstant * pE / (TXC_K_ONE_HALF + pE);
		bP = p->productionConstant * bP / (TXC_K_ONE_HALF + bP);
		*baseProdEffect = bP * p->qty;
		*maxProdEffect = pE * p->qty;
	}else{
		pE = pow(pE,TXL_PROD_EXPONENT);
		pE = p->productionConstant * pE / (TXL_K_ONE_HALF + pE);
		bP = p->productionConstant * bP / (TXL_K_ONE_HALF + bP);
		dE = pow(dE,MESS_DECAY_EXPONENT);
		dE = p->decayConstant * dE / (MESS_DECAY_K_ONE_HALF + dE);
		bD = p->decayConstant * bD / (MESS_DECAY_K_ONE_HALF+ bD);
		*baseProdEffect = bP;
		*maxProdEffect = pE;
		*maxDecayEffect = dE;
		*baseDecayEffect = bD;
	}
}
