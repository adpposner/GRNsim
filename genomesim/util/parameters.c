// parameters.c - Functions for parsing xml for network generation
//a few utility functions are used for simulation setup, but not actually part of simulation
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
#include "../include/parameters.h"
#include "../include/parameters_p.h"
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <string.h>
#include "../include/globals.h"

#define LOCAL_XSTR(s)	LOCAL_STR(s)
#define LOCAL_STR(s)	#s

#define FOUND0	1<<0
#define FOUND1	1<<1
#define FOUND2	1<<2
#define FOUND3	1<<3
#define FOUND4	1<<4
#define FOUND5	1<<5
#define FOUND6	1<<6
#define FOUNDALL5 (FOUND0 | FOUND1 | FOUND2 | FOUND3 |FOUND4)


#define XMLNODEASSERTION(node,namestr,docum)	\
	do {	\
		if (strcmp(node->name,(const xmlChar *) namestr)) {  \
			fprintf(stderr,"node name failure - %s != %s\n",node->name,namestr); \
			xmlFreeDoc(docum); \
			exit(-88); \
		} \
	}while(0);

genesDims globalDims;

int networkGenSeed;


static int getXMLint(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

	xmlChar *key;
	int toRet;
	int found = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);

			toRet = strtol((const char *)key, NULL, 10);

			xmlFree(key);
			found = 1;
		}
		cur = cur->next;
	}
	assert(found);
	return toRet;
}


static numeric_t_rp getXMLnumeric(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

	xmlChar *key;
	numeric_t_rp toRet;
	int found = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);

#ifdef DOUBLE_PRECISION_RP
			toRet = strtod((const char *)key, NULL);
#else
			toRet = strtof((const char *)key, NULL);
#endif
			xmlFree(key);
			found = 1;
		}
		cur = cur->next;
	}
	assert(found);
	return toRet;
}

static unsigned char getXMLchar(xmlDocPtr doc, xmlNodePtr cur, const char * name) {
	xmlChar *key;
	unsigned char toRet;
	int found = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar*) name))) {
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			toRet = key[0];
			xmlFree(key);
			found = 1;
		}
		cur = cur->next;
	}
	assert(found);
	return toRet;
}

static void parseDims(xmlDocPtr doc, xmlNodePtr cur) {
	globalDims.nMess = globalDims.nMicro = 0;
	int foundnMess, foundnMicro;
	foundnMicro = foundnMess = 0;
	cur = cur->xmlChildrenNode;
	//make sure it matches the right name
	while (cur != NULL) {
		if (!foundnMess)
		{globalDims.nMess = getXMLint(doc, cur, "nMess"); foundnMess = 1;}
		if (!foundnMicro)
		{globalDims.nMicro = getXMLint(doc, cur, "nMicro"); foundnMicro = 1;}


		cur = cur->next;
	}
	if (foundnMicro && foundnMess) {
		globalDims.nDNA = globalDims.nMess + globalDims.nMicro;
		globalDims.nProt = globalDims.nMess;

		return;
	} else {
		fprintf(stderr, "cannot parse dims\n");
		exit(-895);
	}

}



static struct ConnectionParametersUnit * parseConnectionParameterStruct(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	struct ConnectionParametersUnit *toRet = malloc(sizeof(*toRet));

	int founddl, founddh, foundeffP, foundByInOut;
	founddl = founddh  = foundeffP = foundByInOut = 0;
	char tmp;

	while (cur != NULL) {
		if (!founddl)
		{toRet->deg_low = getXMLint(doc, cur, "degree_low"); founddl = 1;} 
		if (!founddh)
		{toRet->deg_high = getXMLint(doc, cur, "degree_high"); founddh = 1;}
		if (!foundeffP)
		{toRet->effect_prob = getXMLnumeric(doc, cur, "effectProb"); foundeffP = 1;}
		if (!foundByInOut){
					tmp = getXMLchar(doc, cur, "inOutDegree");
					switch(tmp){
					case 'I':
					case 'i': foundByInOut = toRet->byInOrOut = INDEGREE; break;
					case 'O':
					case 'o': foundByInOut = toRet->byInOrOut = OUTDEGREE; break;
					default: exit(-5412);
					}foundByInOut = 1;}
		cur = cur->next;
	}
	if (founddh && founddl && foundeffP && foundByInOut) {
		return toRet;
	} else {
		fprintf(stderr, "failure parsing connparms\n");
		exit(-1000);
	}


}

static void parseTFConnectionParameters(xmlDocPtr doc, xmlNodePtr cur){
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if (!xmlStrcasecmp(cur->name, (const xmlChar *) "coding"))
			system_parameters.connParms.tfCodingConns = parseConnectionParameterStruct(doc,cur);
		else if (!xmlStrcasecmp(cur->name, (const xmlChar *) "noncoding"))
			system_parameters.connParms.tfNoncodingConns = parseConnectionParameterStruct(doc,cur);
		cur = cur->next;
	}
}

static void parseConnectionParameters(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "miRNAs")))
			system_parameters.connParms.miRConns = parseConnectionParameterStruct(doc, cur);

		if ((!xmlStrcmp(cur->name, (const xmlChar *) "TFs")))
			parseTFConnectionParameters(doc,cur);

		cur = cur->next;
	}
}

static void parseRandomSeeds(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int foundNetGen;
	foundNetGen = 0;
	while (cur != NULL) {
		if (!foundNetGen)
		{networkGenSeed = getXMLint(doc, cur, "networkGenSeed"); foundNetGen = 1;}
		cur = cur->next;
	}
	if (foundNetGen) {
		return;
	} else {
		fprintf(stderr, "failure parsing seeds\n");
		exit(-1000);
	}
}

static void parseMessRates(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {
		if (!(found & FOUND0)) {
			system_parameters.rates.mess.production = getXMLnumeric(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.mess.decay = getXMLnumeric(doc, cur, "decay");
			found += FOUND1;
		}
		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1)) {
		fprintf(stderr, "failure parsing messrates\n");
		exit(-125);
	}
}

static void parseMicroRates(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {
		if (!(found & FOUND0)) {
			system_parameters.rates.micro.production = getXMLnumeric(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.micro.decay = getXMLnumeric(doc, cur, "decay");
			found += FOUND1;
		}
		if (!(found & FOUND2)) {
			system_parameters.rates.micro.affinity = getXMLnumeric(doc, cur, "affinity");
			found += FOUND2;
		}
		if (!(found & FOUND3)) {
			system_parameters.rates.micro.effectStrength = getXMLnumeric(doc, cur, "effectStrength");
			found += FOUND3;
		}
		if (!(found & FOUND4)) {
			system_parameters.rates.micro.frequency = getXMLnumeric(doc, cur, "frequency");
			found += FOUND4;
		}
		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2 | FOUND3 | FOUND4)) {
		fprintf(stderr, "failure parsing microrates\n");
		exit(-124);
	}
}

static void parseTFRates(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {
		if (!(found & FOUND0)) {
			system_parameters.rates.tf.production = getXMLnumeric(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.tf.decay = getXMLnumeric(doc, cur, "decay");
			found += FOUND1;
		}
		if (!(found & FOUND2)) {
			system_parameters.rates.tf.affinity = getXMLnumeric(doc, cur, "affinity");
			found += FOUND2;
		}
		if (!(found & FOUND3)) {
			system_parameters.rates.tf.effectMin = getXMLnumeric(doc, cur, "effectMin");
			found += FOUND3;
		}
		if (!(found & FOUND4)) {
			system_parameters.rates.tf.effectMax = getXMLnumeric(doc, cur, "effectMax");
			found += FOUND4;
		}
		if (!(found & FOUND5)) {
			system_parameters.rates.tf.frequency = getXMLnumeric(doc, cur, "frequency");
			found += FOUND5;
		}
		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2 | FOUND3 | FOUND4 | FOUND5)) {
		fprintf(stderr, "failure parsing tfrates\n");
		exit(-12);
	}
}


static void parseBaseRates(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {
		if (!(found & FOUND0))
			if (!xmlStrcmp(cur->name, (const xmlChar * )"mess")) {
				parseMessRates(doc, cur);
				found += FOUND0;
			}
		if (!(found & FOUND1))
			if (!xmlStrcmp(cur->name, (const xmlChar *) "micro")) {
				parseMicroRates(doc, cur);
				found += FOUND1;
			}
		if (!(found & FOUND2))
			if (!xmlStrcmp(cur->name, (const xmlChar *) "TF")) {
				parseTFRates(doc, cur);
				found += FOUND2;
			}

		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2)) {
		fprintf(stderr, "failure parsing connparms\n");
		exit(-1000);
	}
	return;


}

static void parseInitialQuantities(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {

		if (!(found & FOUND0)) {
			system_parameters.initialQtys.coding = getXMLint(doc, cur, "coding");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.initialQtys.noncoding = getXMLint(doc, cur, "noncoding");
			found += FOUND1;
		}
		if (!(found & FOUND2)) {
			system_parameters.initialQtys.mess = getXMLint(doc, cur, "mess");
			found += FOUND2;
		}
		if (!(found & FOUND3)) {
			system_parameters.initialQtys.micro = getXMLint(doc, cur, "micro");
			found += FOUND3;
		}
		if (!(found & FOUND4)) {
			system_parameters.initialQtys.TF = getXMLint(doc, cur, "TF");
			found += FOUND4;
		}
		cur = cur->next;
	}
	if (found == FOUNDALL5) return;
	else {
		fprintf(stderr, "failure parsing initialqtys\n");
		exit(-19000);
	}
}



#define CONNPARMWARNINGFMTSTRING  "Warning: %s connection degree %s for by %s = %d > %s = %zd\t defaulting to %s\n"

static void checkConnParmsValidity_(struct ConnectionParametersUnit * parms, genesDims *gD) {
	ulong_type hiVal;
	if (parms == NULL)
	{fprintf(stderr,"not initialized parms");exit(-272);}
	const char * hiValName;
	if (parms->connTp == PROTEIN) {
		if (parms->byInOrOut == INDEGREE)
		{hiVal = gD->nProt; hiValName = "nProt";}
		else
		{hiVal = gD->nDNA; hiValName = "nDNA";}
	}
	else if (parms->connTp == MICRO) {
		if (parms->byInOrOut == INDEGREE)
		{hiVal = gD->nMicro; hiValName = "nMicro";}
		else
		{hiVal = gD->nMess; hiValName = "nMess";}
	}
	const char * degType = (parms->byInOrOut == OUTDEGREE) ? "Outdegree" : "Indegree";
	if (parms->deg_high > hiVal){
		fprintf(stderr,CONNPARMWARNINGFMTSTRING,
			parms->connTp == MICRO ? "miRNA" : "TF","high",degType,parms->deg_high,hiValName,hiVal,hiValName);
		exit(-782);
	}
	if (parms->deg_low > hiVal){
		fprintf(stderr,CONNPARMWARNINGFMTSTRING,
			parms->connTp == MICRO ? "miRNA" : "TF","low",degType,parms->deg_low,hiValName,hiVal,hiValName);
		exit(-783);	
	}

}


static void check_parameters_validity() {
	checkConnParmsValidity_(system_parameters.connParms.miRConns, &globalDims);
	checkConnParmsValidity_(system_parameters.connParms.tfCodingConns, &globalDims);
	checkConnParmsValidity_(system_parameters.connParms.tfNoncodingConns,&globalDims);
}


void initializeParameters(const char * cwd) {
	xmlDocPtr doc;
	xmlNodePtr cur;
	char str[400];
	sprintf(str, "%s/config.xml", cwd);
	doc = xmlParseFile(str);
	system_parameters.connParms.miRConns = system_parameters.connParms.tfCodingConns = 
	system_parameters.connParms.tfNoncodingConns = NULL;
	if (doc == NULL) {
		fprintf(stderr, "config.xml parsing failed");
		exit(-85);
	}

	cur = xmlDocGetRootElement(doc);

	if (cur == NULL) {
		fprintf(stderr, "empty document\n");
		xmlFreeDoc(doc);
		exit(-84);
	}

	if (xmlStrcmp(cur->name, (const xmlChar *) "configuration")) {
		fprintf(stderr, "Root node failure - should be \"configuration\", read as \"%s\"\n", cur->name);
		xmlFreeDoc(doc);
		exit(-86);
	}
	cur = cur->xmlChildrenNode;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "globalDimensions")))
			parseDims(doc, cur);
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "connectionParameters")))
			parseConnectionParameters(doc, cur);
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "randomSeeds")))
			parseRandomSeeds(doc, cur);
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "baseRates")))
			parseBaseRates(doc, cur);
		if ((!xmlStrcmp(cur->name, (const xmlChar *) "initialQuantities")))
			parseInitialQuantities(doc, cur);
		cur = cur->next;

	}
	//Check validity
	check_parameters_validity();
	xmlFreeDoc(doc);

}

void destroyParameters() {
	free(system_parameters.connParms.miRConns);
	free(system_parameters.connParms.tfCodingConns);
	free(system_parameters.connParms.tfNoncodingConns);
}


// void printParms() {
// 	printf("globalDims\n\tnMess=%zd\n\tnMicro=%zd\n", globalDims.nMess, globalDims.nMicro);
// 	printf("connectionParameters\n");
// 	printf("\tmiRNAs\n\t\tdeg_low=%d\n\t\tdeg_high=%d\n\t\teffect_prob=%f\n", miRparms.deg_low, miRparms.deg_high, miRparms.effect_prob);
// 	printf("\tTFs\n\t\tdeg_low=%d\n\t\tdeg_high=%d\n\t\teffect_prob=%f\n", TFparms.deg_low, TFparms.deg_high, TFparms.effect_prob);
// 	//printf("system_parameters.rates\n\tproteinProd=%f\n\tdefaultRate=%f\n",system_parameters.rates.proteinProd,system_parameters.rates.defaultRate);
// 	printf("randomSeeds\n");
// 	printf("\tnetworkGenSeed=%d\n\tsimulationSeed=%d\n", networkGenSeed, simulationSeed);
// }




//The following several functions provide default values for assoc/dissoc if none given
effect_t_rp miREffectStrength(effect_t_rp *r){return (r) ? *r : system_parameters.rates.micro.effectStrength;}
rate_t_rp TFAffinity(rate_t_rp *r){return (r) ? *r : system_parameters.rates.tf.affinity;}
#include "../include/randomnumbers.h"
effect_t_rp TFEffectStrength(effect_t_rp *r){if (r) return *r; 
	effect_t_rp effect;
	 randomUniformNumberArray(1, system_parameters.rates.tf.effectMin, 
	 	system_parameters.rates.tf.effectMax, &effect);
	return effect;}

rate_t_rp messengerProductionRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.mess.production;
}
rate_t_rp microProductionRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.micro.production;
}
rate_t_rp proteinProductionRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.tf.production;
}
rate_t_rp messengerDecayRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.mess.decay;
}
rate_t_rp microDecayRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.micro.decay;
}
rate_t_rp proteinDecayRate(rate_t_rp *r){
	return (r) ? *r : system_parameters.rates.tf.decay;
}


rate_t_rp miRAssocConst(rate_t_rp *assoc){
	if (assoc)
		return *assoc;
	rate_t_rp scale = system_parameters.rates.micro.frequency;
	rate_t_rp aff = system_parameters.rates.micro.affinity;
	return (2.0*scale*aff)/(1.0+aff);
}

rate_t_rp miRDissocConst(rate_t_rp * dissoc){
	if (dissoc)
		return *dissoc;
	rate_t_rp scale = system_parameters.rates.micro.frequency;
	rate_t_rp aff = system_parameters.rates.micro.affinity;
	return (2.0*scale)/(1.0+aff);
}

rate_t_rp TFAssocConst(rate_t_rp *assoc){
	if (assoc)
		return *assoc;
	rate_t_rp scale = system_parameters.rates.tf.frequency;
	rate_t_rp aff = system_parameters.rates.tf.affinity;
	return (2.0*scale*aff)/(1.0+aff);
}

rate_t_rp TFDissocConst(rate_t_rp * dissoc){
	if (dissoc)
		return *dissoc;
	rate_t_rp scale = system_parameters.rates.tf.frequency;
	rate_t_rp aff = system_parameters.rates.tf.affinity;
	return (2.0*scale)/(1.0+aff);
}

unsigned char connsByOutdegree(species_t sp){
	if ((sp == MICRO) || (sp == MESSENGER)){
		return system_parameters.connParms.miRConns->byInOrOut;
	}else if (sp & CODING){
		return system_parameters.connParms.tfCodingConns->byInOrOut;
	}else if (sp& NONCODING) {
		return system_parameters.connParms.tfNoncodingConns->byInOrOut;
	}else{
		exit(-2329);
	}
}

UnsignedIntArray getRandomConnectionNumbersForElements(unsigned int leftLen, unsigned int rightLen, species_t left,species_t right) {
	struct ConnectionParametersUnit * parms;

	if (left & CODING){
		assert(right == PROTEIN);
		parms = system_parameters.connParms.tfCodingConns;
	} else if (left & NONCODING){
		assert(right == PROTEIN);
		parms = system_parameters.connParms.tfNoncodingConns;
	}else{
		assert(left == MESSENGER);
		assert(right == MICRO);
		parms = system_parameters.connParms.miRConns;
	}
	UnsignedIntArray toRet;
	unsigned int nElems = (parms->byInOrOut == INDEGREE) ? leftLen : rightLen;
	toRet = UnsignedIntArray_alloc(nElems);
	if ((left == MESSENGER) && BIMODAL_MIRNA_RP){

		randomBimodalBinomialArray(nElems, parms->deg_low, parms->deg_high, parms->effect_prob, 
			BIMODAL_P2_RP * parms->effect_prob, BIMODAL_SELECT_RP, toRet.data);
	}else{

		randomBinomialArray(nElems, parms->deg_low, parms->deg_high, parms->effect_prob, toRet.data);
	}
	return toRet;
}


ulong_type getInitialQty(species_t sp){
	switch(sp){
		case CODING: return system_parameters.initialQtys.coding;
		case NONCODING: return system_parameters.initialQtys.noncoding;
		case MESSENGER: return system_parameters.initialQtys.mess;
		case MICRO: return system_parameters.initialQtys.micro;
		case PROTEIN: return system_parameters.initialQtys.TF;
		case MESSMIR:
		case TFDNA: return 0;
		default: exit(-267398);
	}
}
