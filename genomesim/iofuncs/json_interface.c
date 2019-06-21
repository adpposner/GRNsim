//json_interface.c - source for JSON I/O, using cJSON library
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
#include "../include/models/basicgeneticelement_p.h"
#include "../include/models/basicgeneticelement.h"

#include "../include/models/geneslist_p.h"
//#include "../include/models/reaction.h"
#include "../include/iofuncs/json_interface.h"
//#include "../include/models/geneslist.h"
#include "../include/globals.h"

#include <string.h>

struct cJSON * UnsignedIntArray_toJSON(UnsignedIntArray *arr);
UnsignedIntArray UnsignedIntArray_fromJSON(struct cJSON * obj);

cJSON * UnsignedIntArray_toJSON(UnsignedIntArray *arr){
	cJSON * toRet;
	toRet = cJSON_CreateIntArray((const int *)arr->data, arr->length);
	return toRet;
}

UnsignedIntArray UnsignedIntArray_fromJSON(cJSON * obj) {
 	UnsignedIntArray toRet;
 	assert(cJSON_IsArray(obj));
 	int len;
 	len = cJSON_GetArraySize(obj);

 	toRet = UnsignedIntArray_alloc(len);
 	unsigned * data = toRet.data;
 	cJSON * item;
 	cJSON_ArrayForEach(item, obj) *data++ = item->valueint;
 	return toRet;
 }


#define InfoToJson(elem,jsonObj) cJSON * item = cJSON_CreateObject(); \
	cJSON_AddStringToObject(item, "name", elem->name);					\
	cJSON_AddNumberToObject(item, "species", elem->species);			\
	cJSON_AddNumberToObject(item, "qty",elem->qty);						\
	cJSON_AddNumberToObject(item, "id", elem->id);						\
	cJSON_AddNumberToObject(item, "isEnabled", elem->isEnabled);		\
	cJSON_AddItemToObject(jsonObj, "info", item);

#define InfoFromJson(elem,obj) cJSON * item = cJSON_GetObjectItem(obj, "info");	\
	snprintf(elem.name,MAXGENENAMELENGTH,"%s",cJSON_GetObjectItem(item, "name")->valuestring); \
	elem.id = cJSON_GetObjectItem(item, "id")->valueint;	\
	elem.qty = cJSON_GetObjectItem(item, "qty")->valueint;	\
	elem.species = (species_t) cJSON_GetObjectItem(item, "species")->valueint; \
	elem.isEnabled = cJSON_GetObjectItem(item, "isEnabled")->valueint;




static UnsignedIntArray getBoundIDsForProducer(Producer * p){
	UnsignedIntArray toRet;
	toRet = UnsignedIntArray_alloc(p->boundelts.length);
	BoundElementArrayIters bIt = getBoundElementArrayIters(&p->boundelts);
	UnsignedIntArrayIters uIt = getUnsignedIntArrayIters(&toRet);
	for(bIt.curr=bIt.start,uIt.curr=uIt.start;bIt.curr!=bIt.end;bIt.curr++,uIt.curr++)
		*uIt.curr = bIt.curr->id;
	assert(uIt.curr == uIt.end);
	return toRet;
}

static UnsignedIntArray getBoundIDsForModulator(Modulator * p){
	UnsignedIntArray toRet;
	toRet = UnsignedIntArray_alloc(p->boundelts.length);
	BoundElementPtrArrayIters bpIt = getBoundElementPtrArrayIters(&p->boundelts);
	UnsignedIntArrayIters uIt = getUnsignedIntArrayIters(&toRet);
	for(bpIt.curr=bpIt.start,uIt.curr=uIt.start;bpIt.curr!=bpIt.end;bpIt.curr++,uIt.curr++)
		*uIt.curr = (*bpIt.curr)->id;
	assert(uIt.curr == uIt.end);
	return toRet;
}

static cJSON * producerBounds_toJSON(Producer *p){
	UnsignedIntArray bds = getBoundIDsForProducer(p);
	cJSON * toRet =cJSON_CreateArray();
	UnsignedIntArrayIters uIt = getUnsignedIntArrayIters(&bds);
	ARRAY_TYPE_FOREACH(uIt) cJSON_AddItemToArray(toRet, cJSON_CreateNumber(*uIt.curr));
	UnsignedIntArray_free(&bds);
	return toRet;
}

static cJSON * modulatorBounds_toJSON(Modulator *m){
	UnsignedIntArray bds = getBoundIDsForModulator(m);
	cJSON * toRet =cJSON_CreateArray();
	UnsignedIntArrayIters uIt = getUnsignedIntArrayIters(&bds);
	ARRAY_TYPE_FOREACH(uIt) cJSON_AddItemToArray(toRet, cJSON_CreateNumber(*uIt.curr));
	UnsignedIntArray_free(&bds);
	return toRet;
}


cJSON * producer_toJSON(Producer *elem){
	cJSON * toRet;
	toRet = cJSON_CreateObject();
	InfoToJson(elem, toRet)
	cJSON_AddNumberToObject(toRet, "productionConstant", elem->productionConstant);
	cJSON_AddNumberToObject(toRet, "decayConstant", elem->decayConstant);
	cJSON_AddNumberToObject(toRet, "produces",getProducesID(elem));
	cJSON_AddItemToObject(toRet, "boundelts", producerBounds_toJSON(elem));
	return toRet;
}




static Producer producerBase_fromJSON(cJSON *obj){
	Producer toRet;
	InfoFromJson(toRet, obj);
	toRet.productionConstant = cJSON_GetObjectItem(obj, "productionConstant")->valuedouble;
	toRet.decayConstant = cJSON_GetObjectItem(obj, "decayConstant")->valuedouble;
	toRet.produces = emptyPmbPtr();
	toRet.boundelts = emptyBoundElementArray();
	toRet.occupancies = emptyOccupancyVector();
	toRet.reactions = nullProducerRxns();
	setProducerQty(&toRet, toRet.qty);
	return toRet;	
}

static Modulator modulatorBase_fromJSON(cJSON *obj){
	Modulator toRet;
	InfoFromJson(toRet, obj);
	toRet.decayConstant = cJSON_GetObjectItem(obj, "decayConstant")->valuedouble;
	toRet.boundelts = emptyBoundElementPtrArray();
	toRet.bindings = emptyReactionPtrArray();
	return toRet;
}

static BoundElement boundBase_fromJSON(cJSON *obj){
	BoundElement toRet;
	InfoFromJson(toRet, obj);
	toRet.left = NULL;
	toRet.right = NULL;
	toRet.producerArrayPos = 0;
	toRet.unbinding = NULL;
	toRet.assocConstant = cJSON_GetObjectItem(obj, "assocConstant")->valuedouble;
	toRet.dissocConstant = cJSON_GetObjectItem(obj, "dissocConstant")->valuedouble;
	toRet.effectStrength = cJSON_GetObjectItem(obj, "effectStrength")->valuedouble;
	#ifdef BOUNDDECAY
	toRet.decayConstant = cJSON_GetObjectItem(obj, "decayConstant")->valuedouble;
	#endif
	return toRet;
}

static ProducerArray producerArrayBase_fromJSON(cJSON *arr){
	ProducerArray toRet;
	cJSON *item;
	ulong_type len = cJSON_GetArraySize(arr);
	toRet = ProducerArray_alloc(len);
	ProducerArrayIters pIt = getProducerArrayIters(&toRet);
	pIt.curr=pIt.start;
	ulong_type k;
	for(k=0;k<len;k++){
		item = cJSON_GetArrayItem(arr, k);
		*pIt.curr++ = producerBase_fromJSON(item);
	}
	assert(pIt.curr==pIt.end);
	return toRet;
}

static ModulatorArray modulatorArrayBase_fromJSON(cJSON *arr){
	ModulatorArray toRet;
	cJSON *item;
	ulong_type len = cJSON_GetArraySize(arr);

	toRet = ModulatorArray_alloc(len);
	ModulatorArrayIters mIt = getModulatorArrayIters(&toRet);
	mIt.curr=mIt.start;
	cJSON_ArrayForEach(item, arr){
		*mIt.curr++ = modulatorBase_fromJSON(item);
	}
	assert(mIt.curr==mIt.end);
	return toRet;
}

static BoundElementArray boundArrayBase_fromJSON(cJSON *arr){
	BoundElementArray toRet;
	cJSON *item;
	ulong_type len = cJSON_GetArraySize(arr);

	toRet = BoundElementArray_alloc(len);
	BoundElementArrayIters bIt = getBoundElementArrayIters(&toRet);
	bIt.curr=bIt.start;
	cJSON_ArrayForEach(item, arr){
		*bIt.curr++ = boundBase_fromJSON(item);
	}
	assert(bIt.curr==bIt.end);
	return toRet;
}



static cJSON * modulator_toJSON(Modulator *elem){
	cJSON * toRet;
	toRet = cJSON_CreateObject();
	InfoToJson(elem, toRet);
	cJSON_AddNumberToObject(toRet, "decayConstant", elem->decayConstant);
	cJSON_AddItemToObject(toRet, "boundelts", modulatorBounds_toJSON(elem));
	return toRet;
}


static cJSON * bound_toJSON(BoundElement *bd){
	cJSON * toRet = cJSON_CreateObject();
	InfoToJson(bd, toRet);
	cJSON_AddNumberToObject(toRet, "left", bd->left->id);
	cJSON_AddNumberToObject(toRet, "right", bd->right->id);
	cJSON_AddNumberToObject(toRet, "assocConstant", bd->assocConstant);
	cJSON_AddNumberToObject(toRet, "dissocConstant", bd->dissocConstant);
	cJSON_AddNumberToObject(toRet, "effectStrength", bd->effectStrength);
	#ifdef BOUNDDECAY
	cJSON_AddNumberToObject(toRet, "decayConstant", bd->decayConstant);
	#endif
	return toRet;
}

static cJSON * producerArray_toJSON(ProducerArray * p,SkipDisabled skipDisabled){
	cJSON * toRet = cJSON_CreateArray();
	ProducerArrayIters pIt = getProducerArrayIters(p);
	ARRAY_TYPE_FOREACH(pIt) {
		if(skipDisabled == SKIP_DISABLED)
			if(pIt.curr->isEnabled == 0)
				continue;
		cJSON_AddItemToArray(toRet, producer_toJSON(pIt.curr));
	}
	return toRet;
}

static cJSON * modulatorArray_toJSON(ModulatorArray * m,SkipDisabled skipDisabled){
	cJSON * toRet = cJSON_CreateArray();
	ModulatorArrayIters mIt = getModulatorArrayIters(m);
	ARRAY_TYPE_FOREACH(mIt) {
		if(skipDisabled == SKIP_DISABLED)
			if(mIt.curr->isEnabled == 0)
				continue;
		cJSON_AddItemToArray(toRet, modulator_toJSON(mIt.curr));
	}
	return toRet;
}

static cJSON * boundArray_toJSON(BoundElementArray * b,SkipDisabled skipDisabled){
	cJSON * toRet = cJSON_CreateArray();
	BoundElementArrayIters bIt = getBoundElementArrayIters(b);
	ARRAY_TYPE_FOREACH(bIt)	{
		if(skipDisabled == SKIP_DISABLED)
			if(bIt.curr->isEnabled == 0)
				continue;
		cJSON_AddItemToArray(toRet, bound_toJSON(bIt.curr));
	}
	return toRet;
}


static cJSON * genesListDims_toJSON(GenesList * g,SkipDisabled skipDisabled){
	cJSON * toRet = cJSON_CreateObject();
	cJSON_AddNumberToObject(toRet, "nDNA", (skipDisabled == INCLUDE_DISABLED) ?  g->nDNA : nEnabledForType(g, DNA));
	cJSON_AddNumberToObject(toRet, "nMess",  (skipDisabled == INCLUDE_DISABLED) ?  g->nMess: nEnabledForType(g, MESSENGER));
	cJSON_AddNumberToObject(toRet, "nMicro",  (skipDisabled == INCLUDE_DISABLED) ?  g->nMicro : nEnabledForType(g, MICRO));
	cJSON_AddNumberToObject(toRet, "nProt",  (skipDisabled == INCLUDE_DISABLED) ?  g->nProt : nEnabledForType(g, PROTEIN));
	cJSON_AddNumberToObject(toRet, "nMessMir",  (skipDisabled == INCLUDE_DISABLED) ?  g->nMessMir: nEnabledForType(g, MESSMIR));
	cJSON_AddNumberToObject(toRet, "nTFDNA",  (skipDisabled == INCLUDE_DISABLED) ?  g->nTFDNA : nEnabledForType(g, TFDNA));
	return toRet;
}

static GenesList *genesListInit_fromJSON(cJSON * obj){
	cJSON * dims = cJSON_GetObjectItem(obj, "dims");
	ulong_type nDNA,nMess,nMicro,nProt,nMessMir,nTFDNA;
	nDNA = nMess = nMicro = nProt = nMessMir = nTFDNA = 0;
	nDNA = cJSON_GetObjectItem(dims, "nDNA")->valueint;
	nMess = cJSON_GetObjectItem(dims, "nMess")->valueint;
	nMicro = cJSON_GetObjectItem(dims, "nMicro")->valueint;
	nProt = cJSON_GetObjectItem(dims, "nProt")->valueint;
	nMessMir = cJSON_GetObjectItem(dims, "nMessMir")->valueint;
	nTFDNA = cJSON_GetObjectItem(dims, "nTFDNA")->valueint;
	assert(nDNA == (nMess + nMicro));
	assert(nProt == nMess);
	GenesList *toRet = genesList_alloc(nMess, nMicro, nMessMir, nTFDNA,0);
	return toRet;
}


static cJSON * genesList_toJSON(GenesList *g,SkipDisabled skipDisabled) {
	cJSON * toRet = cJSON_CreateObject();
	cJSON_AddItemToObject(toRet, "dims", genesListDims_toJSON(g,skipDisabled));
	cJSON_AddItemToObject(toRet, "producers", producerArray_toJSON(&g->producers,skipDisabled));
	cJSON_AddItemToObject(toRet, "modulators", modulatorArray_toJSON(&g->modulators,skipDisabled));
	cJSON_AddItemToObject(toRet, "boundelts", boundArray_toJSON(&g->bounds,skipDisabled));
	return toRet;
}







ulong_type getMaxIndex(GenesList *g){
	ulong_type toRet = 0;
	ProducerArrayIters pIt = getProducerArrayIters(&g->producers);
	ModulatorArrayIters mIt = getModulatorArrayIters(&g->modulators);
	BoundElementArrayIters bIt = getBoundElementArrayIters(&g->bounds);
	ARRAY_TYPE_FOREACH(pIt)	if(pIt.curr->id > toRet) toRet = pIt.curr->id;
	ARRAY_TYPE_FOREACH(mIt)	if(mIt.curr->id > toRet) toRet = mIt.curr->id;
	ARRAY_TYPE_FOREACH(bIt)	if(bIt.curr->id > toRet) toRet = bIt.curr->id;
	return toRet;
}



pmbPtrArray buildIDIndex(GenesList *g){
	ulong_type maxIndex = getMaxIndex(g);
	pmbPtrArray toRet = pmbPtrArray_alloc(maxIndex+1);
	GenesListIterator gIt = getGenesListIters(g);
	ARRAY_TYPE_FOREACH(gIt.pIt)	{
		toRet.data[gIt.pIt.curr->id] = (pmbPtr) {.producer = gIt.pIt.curr,
			.species = gIt.pIt.curr->species};}
	ARRAY_TYPE_FOREACH(gIt.mIt)	{
		toRet.data[gIt.mIt.curr->id] = (pmbPtr) {.modulator = gIt.mIt.curr,
			.species = gIt.mIt.curr->species};}
	ARRAY_TYPE_FOREACH(gIt.bIt)	{
		toRet.data[gIt.bIt.curr->id] = (pmbPtr) {.bound = gIt.bIt.curr,
			.species = gIt.bIt.curr->species};}
	return toRet;	
}


static void connectProduces(Producer *p,cJSON * pObj,pmbPtrArray *idx) {
	//first connect
	cJSON * info = cJSON_GetObjectItem(pObj, "info");
	ulong_type idObj = cJSON_GetObjectItem(info, "id")->valueint;
	assert(p->id == idObj);
	ulong_type producesID = cJSON_GetObjectItem(pObj, "produces")->valueint;
	p->produces = idx->data[producesID];
}

static void connectLeftRight(BoundElement *b, cJSON * bObj, pmbPtrArray *idx) {
	ulong_type lId, rId;
	cJSON * info = cJSON_GetObjectItem(bObj, "info");
	ulong_type id = cJSON_GetObjectItem(info, "id")->valueint;
	assert(b->id == id);
	lId = cJSON_GetObjectItem(bObj, "left")->valueint;
	rId = cJSON_GetObjectItem(bObj, "right")->valueint;
	b->left = idx->data[lId].producer;
	b->right = idx->data[rId].modulator;
}


static GenesList * genesListBase_fromJSON(cJSON * jdata){

	GenesList * toRet = genesListInit_fromJSON(jdata);
	cJSON * pArr = cJSON_GetObjectItem(jdata, "producers");
	toRet->producers = producerArrayBase_fromJSON(pArr);
	cJSON * mArr = cJSON_GetObjectItem(jdata, "modulators");
	toRet->modulators = modulatorArrayBase_fromJSON(mArr);
	cJSON * bArr = cJSON_GetObjectItem(jdata, "boundelts");
	toRet->bounds = boundArrayBase_fromJSON(bArr);
	pmbPtrArray idx = buildIDIndex(toRet);

	GenesListIterator g = getGenesListIters(toRet);
	ulong_type k =0;
	cJSON *pObj,*bObj;
	ARRAY_TYPE_FOREACH(g.pIt){
		pObj = cJSON_GetArrayItem(pArr, k++);
		connectProduces(g.pIt.curr, pObj, &idx);
	}
	assert(k==toRet->producers.length);
	k=0;
	ARRAY_TYPE_FOREACH(g.bIt){
		bObj = cJSON_GetArrayItem(bArr, k++);
		connectLeftRight(g.bIt.curr, bObj, &idx);
	}
	assert(k==toRet->bounds.length);
	assignBoundElements(toRet);
	//cleanup
	pmbPtrArray_free(&idx);
	return toRet;

}

void writeGenesList_JSON(GenesList * g,const char * filename, SkipDisabled skipDisabled){
	cJSON * jsonout;
	char * strOut;
	jsonout = genesList_toJSON(g,skipDisabled);
	strOut = cJSON_Print(jsonout);
	FILE * outFile;
	outFile = fopen(filename,"wb");
	printf("%s\n",filename);
	fwrite(strOut,sizeof(*strOut),strlen(strOut),outFile);
	fclose(outFile);
	cJSON_Delete(jsonout);
	free(strOut);
}

GenesList * readGenesList_JSON(const char * filename){
	cJSON * jsondata;
	ulong_type length;
	char * strIn;
	FILE * fp = fopen(filename,"rb");
	fseek(fp, 0, SEEK_END);
	length = ftell(fp);
	fseek(fp,0,SEEK_SET);
	strIn = malloc(length);
	fread(strIn,1,length,fp);
	fclose(fp);
	jsondata = cJSON_Parse(strIn);
	free(strIn);
	GenesList *toRet = genesListBase_fromJSON(jsondata);
	extractInitialQuantities(toRet);
	cJSON_Delete(jsondata);
	return toRet;
}

