#include "../include/globals.h"
#include "../include/parameters.h"
#include "../include/models/geneslist.h"
#include "../include/models/geneslist_p.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//fix this
#include "../include/models/basicgeneticelement.h"
#include "../include/models/basicgeneticelement_p.h"



GenesListIterator getGenesListIters(GenesList * g){
	ProducerArrayIters pIt = getProducerArrayIters(&g->producers);
	ModulatorArrayIters mIt = getModulatorArrayIters(&g->modulators);
	BoundElementArrayIters bIt = getBoundElementArrayIters(&g->bounds);
	GenesListIterator toRet = {.pIt = pIt, .mIt = mIt, .bIt = bIt};
	toRet.pIt.curr=toRet.pIt.start;
	toRet.mIt.curr=toRet.mIt.start;
	toRet.bIt.curr=toRet.bIt.start;
	return toRet;
}

ulong_type nMessengers(GenesList * g){return g->nMess;}

ulong_type nMicro_total(GenesList * g){ return g->nMicro;}


ulong_type numElemsForType(GenesList *g,species_t sptype,SkipDisabled skipDisabled){
	ulong_type toRet = 0;
	GenesListIterator gIt = getGenesListIters(g);
	ARRAY_TYPE_FOREACH(gIt.pIt) {
		if (skipDisabled == SKIP_DISABLED)
			if (gIt.pIt.curr->isEnabled == 0) continue;
		if(gIt.pIt.curr->species & sptype) toRet++;
	}
	ARRAY_TYPE_FOREACH(gIt.mIt)	{
		if (skipDisabled == SKIP_DISABLED)
			if (gIt.mIt.curr->isEnabled == 0) continue;
		if(gIt.mIt.curr->species & sptype) toRet++;}
	ARRAY_TYPE_FOREACH(gIt.bIt) {
		if (skipDisabled == SKIP_DISABLED)
			if (gIt.bIt.curr->isEnabled == 0) continue;
		if(gIt.bIt.curr->species & sptype) toRet++;}
	//Should know the number if we are including disabled, this
	//is just a check
	if (skipDisabled == INCLUDE_DISABLED)
		switch(sptype){
			case MESSENGER:
			case PROTEIN:
			case CODING: assert(toRet == g->nMess);break;
			case MICRO:
			case NONCODING: assert(toRet == g->nMicro);break;
			case DNA: assert(toRet == g->nDNA);break;
			case MESSMIR: assert(toRet == g->nMessMir);break;
			case TFDNA: assert(toRet == g->nTFDNA);break;
			case PRODUCER_SPECIES: assert(toRet == (g->nDNA + g->nMess));break;
			case MODULATOR_SPECIES: assert(toRet == (g->nMicro + g->nProt));break;
			case BOUND_SPECIES: assert(toRet == (g->nTFDNA + g->nMessMir));break;
			default: perror("Incorrect num elems for type"); exit(-1103);
		}
	return toRet;
}

ulong_type nEnabledForType(GenesList * g,species_t species){
	ulong_type nEnabled = numElemsForType(g, species, SKIP_DISABLED);
	return nEnabled;
}

pmbPtrArray getPtrArrayForType(GenesList *g, species_t sptype){
	pmbPtrArray toRet;
	ulong_type arrLen = numElemsForType(g, sptype,INCLUDE_DISABLED);
	toRet = pmbPtrArray_alloc(arrLen);
	GenesListIterator gIt =  getGenesListIters(g);
	pmbPtrArrayIters pmbIt = getpmbPtrArrayIters(&toRet);
	pmbIt.curr=pmbIt.start;
	ARRAY_TYPE_FOREACH(gIt.pIt)
		if(gIt.pIt.curr->species & sptype)
			*pmbIt.curr++ = (pmbPtr){.producer = gIt.pIt.curr,.species = sptype};
	ARRAY_TYPE_FOREACH(gIt.mIt)
		if (gIt.mIt.curr->species & sptype)
			*pmbIt.curr++ = (pmbPtr){.modulator = gIt.mIt.curr,.species = sptype};
	ARRAY_TYPE_FOREACH(gIt.bIt)
		if (gIt.bIt.curr->species & sptype)
			*pmbIt.curr++ = (pmbPtr){.bound = gIt.bIt.curr,.species = sptype};
	assert(pmbIt.curr==pmbIt.end);
	return toRet;
}



unsigned char isProducerEnabled(Producer * producer){
	if(producer->species == NONCODING)
		 return producer->produces.modulator->isEnabled;
	return 1;
}

unsigned char isModulatorEnabled(Modulator * modulator){
	return modulator->isEnabled;
}

unsigned char isBoundEnabled(BoundElement * bd){
	if(bd->species == MESSMIR)
		return bd->right->isEnabled;
	else if (bd->left->species == NONCODING)
		return bd->left->produces.modulator->isEnabled;
	else return 1;
}

unsigned isEnabled(const pmbPtr pmb) {
	switch (pmb.species) {
	case CODING:
	case PROTEIN:
	case MESSENGER: return 1;
	case TFDNA:
	case MESSMIR: return isBoundEnabled(pmb.bound);
	case MICRO: return isModulatorEnabled(pmb.modulator);
	case NONCODING: return isProducerEnabled(pmb.producer);
	default: perror("Invalid ptr for enabled check"); exit(-2312);
	}
}

void toggleFamilyForDisabled(GenesList *g){
	GenesListIterator it = getGenesListIters(g);
	ARRAY_TYPE_FOREACH(it.pIt) it.pIt.curr->isEnabled = isProducerEnabled(it.pIt.curr);
	ARRAY_TYPE_FOREACH(it.mIt) it.mIt.curr->isEnabled = isModulatorEnabled(it.mIt.curr);
	ARRAY_TYPE_FOREACH(it.bIt) it.bIt.curr->isEnabled = isBoundEnabled(it.bIt.curr);
}



GenesList * genesList_alloc(const ulong_type nMess, const ulong_type nMicro, 
	const ulong_type nMessMir, const ulong_type nTFDNA,unsigned char allocInternalArrays){
	GenesList * toRet = malloc(sizeof(*toRet));
	toRet->nMess = nMess;
	toRet->nMicro = nMicro;
	toRet->nProt = nMess;
	toRet->nDNA = nMess + nMicro;
	toRet->nMessMir = nMessMir;
	toRet->nTFDNA = nTFDNA;
	if(allocInternalArrays){
	toRet->producers = ProducerArray_alloc(toRet->nDNA + toRet->nMess);
	toRet->modulators = ModulatorArray_alloc(toRet->nMicro + toRet->nProt);
	toRet->bounds = BoundElementArray_alloc(toRet->nMessMir + toRet->nTFDNA);
	}else {toRet->producers=emptyProducerArray();toRet->modulators=emptyModulatorArray(); toRet->bounds=emptyBoundElementArray();}
	toRet->ICs = NULL;
	return toRet;
}

void initialQuantities_alloc(GenesList *g){
	g->ICs = malloc(sizeof(*g->ICs));
	g->ICs->prodICs = UnsignedIntArray_alloc(g->producers.length);
	g->ICs->modICs = UnsignedIntArray_alloc(g->modulators.length);
	g->ICs->boundICs = UnsignedIntArray_alloc(g->bounds.length);
	memset(g->ICs->prodICs.data, 0, g->ICs->prodICs.length * sizeof(*g->ICs->prodICs.data));
	memset(g->ICs->modICs.data, 0, g->ICs->modICs.length * sizeof(*g->ICs->modICs.data));
	memset(g->ICs->boundICs.data, 0, g->ICs->boundICs.length * sizeof(*g->ICs->boundICs.data));
}

void initialQuantities_free(initialQuantities * ic){
	if (ic)
	{
		UnsignedIntArray_free(&ic->prodICs);
		UnsignedIntArray_free(&ic->modICs);
		UnsignedIntArray_free(&ic->boundICs);
	}
}

void genesList_free(GenesList * tofree)
{ 
	ProducerArray_free(&tofree->producers);
	ModulatorArray_free(&tofree->modulators);
	BoundElementArray_free(&tofree->bounds);
	//GenesList emptyList;
	initialQuantities_free(tofree->ICs);
	free(tofree->ICs);
	free(tofree);

}



// static GenesList genesList_base_init(const ulong_type nMess,const ulong_type nMicro)
// {
// 	//GenesList toRet = {.nMess = nMess, .nMicro = nMicro, .nProt = nMess, .nDNA = nMess + nMicro};
// 	GenesList * toRet = genesList_alloc(nMess, nMicro, 0, 0);
// 	return *toRet;
// }

void genesList_bound_quantities_init(GenesList *g, const ulong_type nMessMir, const ulong_type nTFDNA){
	assert(g->nMessMir == 0);
	assert(g->nTFDNA == 0);
	*(ulong_type*) &g->nMessMir = nMessMir;
	*(ulong_type*) &g->nTFDNA = nTFDNA;
	if(!g->bounds.length)
		BoundElementArray_alloc(g->nMessMir + g->nTFDNA);
}



static void assignBoundsToProducers(const ProducerArray *producers, BoundElementArray *bd){
	ProducerArrayIters pIt = getProducerArrayIters(producers);
	sortBoundElements(*bd);
	for(pIt.curr=pIt.start;pIt.curr!=pIt.end;pIt.curr++)
		assignBoundElementsToProducer(pIt.curr, bd);
}

static void assignBoundsToModulators(const ModulatorArray *m, BoundElementArray *bd){
	ModulatorArrayIters mIt = getModulatorArrayIters(m);
	for(mIt.curr=mIt.start;mIt.curr!=mIt.end;mIt.curr++)
		initBoundElementPtrsForModulator(mIt.curr, bd);
}

void assignBoundElements(GenesList *g){
	assignBoundsToProducers(&g->producers, &g->bounds);
	assignBoundsToModulators(&g->modulators, &g->bounds);
}


void initDefaultGenesListQuantities(GenesList *g) {
	GenesListIterator gIt = getGenesListIters(g);
	ARRAY_TYPE_FOREACH(gIt.pIt)
			setDefaultProducerQuantity(gIt.pIt.curr);
	ARRAY_TYPE_FOREACH(gIt.mIt)
			setDefaultModulatorQuantity(gIt.mIt.curr);
	ARRAY_TYPE_FOREACH(gIt.bIt)
			setDefaultBoundQuantity(gIt.bIt.curr);
}


GenesList * genesList_base_generate(ulong_type nMess, ulong_type nMicro)
{
	GenesList *toRet = genesList_alloc(nMess, nMicro, 0, 0,1);
	GenesListIterator gIt = getGenesListIters(toRet);
	ulong_type k;
	for (k = 0; k < nMess; k++){
		initBasicTF(NULL, gIt.pIt.curr,gIt.pIt.curr+1,gIt.mIt.curr++,1);
		gIt.pIt.curr+=2;
	}
	for (k = 0; k < nMicro; k++)
		initBasicMicro(NULL, gIt.pIt.curr++, gIt.mIt.curr++,1);
	return toRet;
}

//Generation of initial quantities presumably from json file
void extractInitialQuantities(GenesList * g){
	if(!g->ICs)
		initialQuantities_alloc(g);
	initialQuantities icTemp = *g->ICs;
	GenesListIterator gIt = getGenesListIters(g);

	UnsignedIntArrayIters picit = getUnsignedIntArrayIters(&icTemp.prodICs);
	UnsignedIntArrayIters modicit = getUnsignedIntArrayIters(&icTemp.modICs);
	UnsignedIntArrayIters bdicit = getUnsignedIntArrayIters(&icTemp.boundICs);
	picit.curr = picit.start;
	modicit.curr=modicit.start;
	bdicit.curr=bdicit.start;
	ARRAY_TYPE_FOREACH(gIt.pIt) *picit.curr++ = gIt.pIt.curr->qty;
	ARRAY_TYPE_FOREACH(gIt.mIt) *modicit.curr++ = gIt.mIt.curr->qty;
	ARRAY_TYPE_FOREACH(gIt.bIt) *bdicit.curr++ = gIt.bIt.curr->qty;
	*(g->ICs)= icTemp;
}

void resetInitialQuantities(GenesList *g){
	GenesListIterator gIt = getGenesListIters(g);
	UnsignedIntArrayIters picit = getUnsignedIntArrayIters(&g->ICs->prodICs);
	UnsignedIntArrayIters modicit = getUnsignedIntArrayIters(&g->ICs->modICs);
	UnsignedIntArrayIters bdicit = getUnsignedIntArrayIters(&g->ICs->boundICs);
	picit.curr = picit.start;
	modicit.curr=modicit.start;
	bdicit.curr=bdicit.start;
	ARRAY_TYPE_FOREACH(gIt.pIt) setProducerQty(gIt.pIt.curr, *picit.curr++);
	ARRAY_TYPE_FOREACH(gIt.mIt) gIt.mIt.curr->qty = *modicit.curr++;
	ARRAY_TYPE_FOREACH(gIt.bIt) gIt.bIt.curr->qty = *bdicit.curr++;
}


