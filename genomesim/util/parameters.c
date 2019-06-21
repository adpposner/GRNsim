#include "../include/parameters.h"
#include "../include/parameters_p.h"
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <string.h>
#include "../include/globals.h"
#include "../include/netgen/networkdistributions.h"
#include "../include/randomnumbers.h"

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

//struct operatingCosts prices;

/*
Cost model parameters inserted here for now


*/

static struct ParameterStruct system_parameters = {.connParms=(connectionParameters){.miRConns = NULL, 
.tfCodingConns = NULL,.tfNoncodingConns=NULL},
.rates = (struct defaultRates){.micro = (struct MicroRates){.production=NULL,
.decay=NULL,.affinity=NULL,.frequency=NULL,.prodEffectStrength=NULL,.decayEffectStrength=NULL},
.mess=(struct MessRates){.production=NULL,.decay=NULL},.tf=(struct TFRates){.production=NULL,
.decay=NULL,.affinity=NULL,.frequency=NULL,.effectSize=NULL},.defaultRate=NULL},
.initialQtys = (struct initialQuantities){.coding=0,.noncoding=0,.mess=0,.micro=0,.TF=0}};



static int getXMLint(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

	xmlChar *key;
	int toRet;
	int found = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            if(!key){
                fprintf(stderr,"%s:%d - failed to get int value for %s\n",__FILE__,__LINE__,name);
                exit(-8754);
            }
			//printf("%s\n",key);
			toRet = strtol((const char *)key, NULL, 10);
			//printf("toRet=%d\n",toRet);
			xmlFree(key);
			found = 1;
		}
		cur = cur->next;
	}
	assert(found);
	return toRet;
}


numeric_t_rp getValueFromSysParam(ConstOrCtsDistribution *unit){
	if (unit->isConst == 1){
        return unit->rate;
    }else if (unit->isConst == -1){
        numeric_t_rp toRet;
        getRandomContinuous(getStream(),1,unit->distn,&toRet);
        return toRet;
    }else{
        fprintf(stderr,"isconst not set for unit\n");
        abort();
    }
}

static numeric_t_rp getXMLnumeric(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

	xmlChar *key;
	numeric_t_rp toRet=0.0;
	int found = 0;
	while (cur != NULL) {
		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
            if(!key){
                fprintf(stderr,"%s:%d - failed to get int value for %s\n",__FILE__,__LINE__,name);
                exit(-8754);
            }
			//printf("%s\n",key);

			toRet = strtod((const char *)key, NULL);

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



static xmlNodePtr getChildByName(xmlDocPtr doc, xmlNodePtr cur, const char * name){
    while(cur!=NULL){
        if(!xmlStrcmp(cur->name,(const xmlChar *)name)){
            return cur;
        }
        cur=cur->next;
    }
    return NULL;
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
	struct ConnectionParametersUnit *toRet = calloc(1,sizeof(*toRet));
    toRet->connectionDistribution = NULL;
    if(!toRet)
    printf("%s:%d malloc fail\n",__FILE__,__LINE__);
	int foundDistribution, foundByInOut;
	foundDistribution = foundByInOut = 0;
	char tmp;

	while (cur != NULL) {
		if (!foundDistribution)
		{
             toRet->connectionDistribution = parseXMLCompoundDiscreteDistribution(doc, cur);
        foundDistribution = 1;
        }
        //printf("founddl%zu\n",toRet.deg_low);}
		
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
	if (foundDistribution && foundByInOut) {
        
		return toRet;
	} else {
		fprintf(stderr, "failure parsing connparms\n");
		exit(-1000);
	}
}

static ConstOrCtsDistribution *getXMLConstOrDist(xmlDocPtr doc, xmlNodePtr cur, const char * name){
            xmlNodePtr effectNode = getChildByName(doc,cur,name);
            
            xmlNodePtr distNode = NULL;
            ConstOrCtsDistribution * toRet;
            toRet = calloc(1,sizeof(*toRet));
            toRet->isConst = 0;

            if (effectNode){
                distNode = getChildByName(doc,effectNode->children,"continuous_distribution");
            if(distNode){
                //is cts distribution
                toRet->isConst = -1;
                toRet->distn= parseXMLCompoundContinuousDistribution(doc,distNode,name);
                return toRet;
            }else{
                 toRet->rate = getXMLnumeric(doc,effectNode,name);
                 toRet->isConst = 1;
                 return toRet;
                }
            }else{
                fprintf(stderr,"Could not find value node for %s\n",name);
                exit(-99);
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
			system_parameters.rates.mess.production = getXMLConstOrDist(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.mess.decay = getXMLConstOrDist(doc, cur, "decay");
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
			system_parameters.rates.micro.production = getXMLConstOrDist(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.micro.decay = getXMLConstOrDist(doc, cur, "decay");
			found += FOUND1;
		}
		if (!(found & FOUND2)) {
			system_parameters.rates.micro.affinity = getXMLConstOrDist(doc, cur, "affinity");
			found += FOUND2;
		}
		if (!(found & FOUND3)) {
			system_parameters.rates.micro.prodEffectStrength = getXMLConstOrDist(doc, cur, "prodEffectStrength");
			found += FOUND3;
		}
		if (!(found & FOUND4)) {
			system_parameters.rates.micro.frequency = getXMLConstOrDist(doc, cur, "frequency");
			found += FOUND4;
		}

        if (!(found & FOUND5)) {
			system_parameters.rates.micro.decayEffectStrength = getXMLConstOrDist(doc, cur, "decayEffectStrength");
			found += FOUND5;
		}		
        cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2 | FOUND3 | FOUND4 | FOUND5)) {
		fprintf(stderr, "failure parsing microrates\n");
		exit(-124);
	}
}

static void parseTFRates(xmlDocPtr doc, xmlNodePtr cur) {
	cur = cur->xmlChildrenNode;
	int found = 0;
	while (cur != NULL) {
		if (!(found & FOUND0)) {
            system_parameters.rates.tf.production = getXMLConstOrDist(doc,cur,"production");
			system_parameters.rates.tf.production = getXMLConstOrDist(doc, cur, "production");
			found += FOUND0;
		}
		if (!(found & FOUND1)) {
			system_parameters.rates.tf.decay = getXMLConstOrDist(doc, cur, "decay");
			found += FOUND1;
		}
		if (!(found & FOUND2)) {
			system_parameters.rates.tf.affinity = getXMLConstOrDist(doc, cur, "affinity");
			found += FOUND2;
		}
        if (!(found & FOUND3)){
            system_parameters.rates.tf.effectSize = getXMLConstOrDist(doc,cur,"effect");
            found += FOUND3;
        }
		
		if (!(found & FOUND4)) {
			system_parameters.rates.tf.frequency = getXMLConstOrDist(doc, cur, "frequency");
			found += FOUND4;
		}
		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2 | FOUND3 | FOUND4)) {
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
        if (!(found & FOUND3))
            if (!xmlStrcmp(cur->name,(const xmlChar *) "defaultRate")){
                system_parameters.rates.defaultRate = getXMLConstOrDist(doc,cur,"defaultRate");
                found += FOUND3;
            }


		cur = cur->next;
	}
	if (found != (FOUND0 | FOUND1 | FOUND2 | FOUND3)) {
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


// static void parsePrices(xmlDocPtr doc, xmlNodePtr cur){
// 	cur=cur->xmlChildrenNode;
// 	int foundPCost,foundMCost,foundMiCost;
// 	foundPCost = foundMCost = foundMiCost = 0;
// 	while (cur !=NULL) {
// 			if(!foundPCost)
// 				{prices.proteinCost = getXMLnumeric(doc,cur,"proteinCost");foundPCost=1;}
// 			if(!foundMCost)
// 				{prices.mRNAcost = getXMLnumeric(doc,cur,"mRNACost");foundMCost=1;}
// 			if(!foundMiCost)
// 				{prices.miRNAcost = getXMLnumeric(doc, cur, "miRNACost");foundMiCost=1;}

// 		cur=cur->next;
// 	}

// 	if (foundPCost && foundMCost && foundMiCost){
// 		return;
// 	}else {
// 		fprintf(stderr,"failure parsing connparms\n");
// 		exit(-1000);
// 	}

// }



// static void parseSystem(xmlDocPtr doc,xmlNodePtr cur) {
// 	cur = cur->xmlChildrenNode;
// 	int foundncycles = 0;
// 	while (cur !=NULL) {
// 		if(!foundmiRsActive)
// 			{miRNAs_active = getXMLint(doc,cur,"miRNAs_active");foundmiRsActive=1;}
// 		if(!foundncycles)
// 			{n_cycles = getXMLint(doc, cur, "n_cycles");foundncycles=1;}
// 		cur = cur->next;
// 	}

// }


static void print_element_names(xmlNode * a_node) {
	xmlNodePtr cur_node;

	for (cur_node = a_node; cur_node; cur_node = cur_node->next) {
		if (cur_node->type == XML_ELEMENT_NODE) {
			//printf("node type: Element,name: %s\n",cur_node->name);
			//key = xmlNodeListGetString(doc,cur_node->xmlChildrenNode,1);
		}
		print_element_names(cur_node->children);
	}
}

#define CONNPARMWARNINGFMTSTRING  "Warning: %s connection degree %s for by %s = %d > %s = %zd\t defaulting to %s\n"

// static void checkConnParmsValidity_(struct ConnectionParametersUnit * parms, genesDims *gD) {
// 	ulong_type hiVal;
// 	if (parms == NULL)
// 	{fprintf(stderr,"not initialized parms");exit(-272);}
// 	const char * hiValName;
// 	if (parms->connTp == PROTEIN) {
// 		if (parms->byInOrOut == INDEGREE)
// 		{hiVal = gD->nProt; hiValName = "nProt";}
// 		else
// 		{hiVal = gD->nDNA; hiValName = "nDNA";}
// 	}
// 	else if (parms->connTp == MICRO) {
// 		if (parms->byInOrOut == INDEGREE)
// 		{hiVal = gD->nMicro; hiValName = "nMicro";}
// 		else
// 		{hiVal = gD->nMess; hiValName = "nMess";}
// 	}
// 	const char * degType = (parms->byInOrOut == OUTDEGREE) ? "Outdegree" : "Indegree";
// 	if (parms->deg_high > hiVal){
// 		fprintf(stderr,CONNPARMWARNINGFMTSTRING,
// 			parms->connTp == MICRO ? "miRNA" : "TF","high",degType,parms->deg_high,hiValName,hiVal,hiValName);
// 		exit(-782);
// 	}
// 	if (parms->deg_low > hiVal){
// 		fprintf(stderr,CONNPARMWARNINGFMTSTRING,
// 			parms->connTp == MICRO ? "miRNA" : "TF","low",degType,parms->deg_low,hiValName,hiVal,hiValName);
// 		exit(-783);	
// 	}

// }


static void check_parameters_validity() {
//	checkConnParmsValidity_(system_parameters.connParms.miRConns, &globalDims);
//	checkConnParmsValidity_(system_parameters.connParms.tfCodingConns, &globalDims);
//	checkConnParmsValidity_(system_parameters.connParms.tfNoncodingConns,&globalDims);
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



// void printParms() {
// 	printf("globalDims\n\tnMess=%zd\n\tnMicro=%zd\n", globalDims.nMess, globalDims.nMicro);
// 	printf("connectionParameters\n");
// 	printf("\tmiRNAs\n\t\tdeg_low=%d\n\t\tdeg_high=%d\n\t\teffect_prob=%f\n", miRparms.deg_low, miRparms.deg_high, miRparms.effect_prob);
// 	printf("\tTFs\n\t\tdeg_low=%d\n\t\tdeg_high=%d\n\t\teffect_prob=%f\n", TFparms.deg_low, TFparms.deg_high, TFparms.effect_prob);
// 	//printf("system_parameters.rates\n\tproteinProd=%f\n\tdefaultRate=%f\n",system_parameters.rates.proteinProd,system_parameters.rates.defaultRate);
// 	printf("randomSeeds\n");
// 	printf("\tnetworkGenSeed=%d\n\tsimulationSeed=%d\n", networkGenSeed, simulationSeed);
// }





#include <openssl/md5.h>



unsigned char * networkParamHash(unsigned char includeNetGenSeed) {
	static unsigned char humanReadableHash[2 * MD5_DIGEST_LENGTH + 1];
	static unsigned char parameterHash[MD5_DIGEST_LENGTH];
	MD5_CTX context;
	char filename[MAX_JSON_FILENAME];
	sprintf(filename,"%s/config.xml",".");
	FILE *configfile= fopen(filename, "r");
	ulong_type length;
	char * strIn;
	FILE * fp = fopen(filename,"rb");
	fseek(fp, 0, SEEK_END);
	length = ftell(fp);
	fseek(fp,0,SEEK_SET);
	strIn = malloc(length+1);
	fread(strIn,1,length,fp);
    strIn[length]='\0';
	fclose(fp);
	
	MD5_Init(&context);
	const char * tags[] = {"globalDimensions","connectionParameters","baseRates","initialQuantities","kinetics"};
	char *pTagStart,*pTagEnd;
	int i;
	for(i=0;i<5;i++){
		pTagStart =strstr(strIn,tags[i]);
		pTagEnd = strstr(pTagStart+1,tags[i]);

		MD5_Update(&context, pTagStart, pTagEnd-pTagStart);
	}
	if (includeNetGenSeed){
		pTagStart=strstr(strIn,"networkGenSeed");
		pTagEnd=strstr(pTagStart+1,"networkGenSeed");
		MD5_Update(&context, pTagStart,pTagEnd-pTagStart);
	}
	
	MD5_Final(parameterHash, &context);
	free(strIn);
	unsigned char letters[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
	for (i = 0; i < MD5_DIGEST_LENGTH - 1; i++)
	{
		humanReadableHash[2 * i] = letters[parameterHash[i] / 52];
		humanReadableHash[2 * i + 1] = letters[parameterHash[i] % 52];
	}
	humanReadableHash[2 * MD5_DIGEST_LENGTH] = 0;
	return humanReadableHash;
}

effect_t_rp miRProdEffectStrength(const effect_t_rp * r){return (r) ? *r : 
getValueFromSysParam(system_parameters.rates.micro.prodEffectStrength);}
effect_t_rp miRDecayEffectStrength(const effect_t_rp * r){return (r) ? *r : 
getValueFromSysParam(system_parameters.rates.micro.decayEffectStrength);}
//effect_t_rp miRProdEffectStrength(effect_t_rp *r){return (r) ? *r : system_parameters.rates.micro.prodEffectStrength;}
//effect_t_rp miRDecayEffectStrength(effect_t_rp *r){return (r) ? *r : system_parameters.rates.micro.decayEffectStrength;}
rate_t_rp TFAffinity(const rate_t_rp *r){return (r) ? *r : getValueFromSysParam(system_parameters.rates.tf.affinity);}
#include "../include/randomnumbers.h"
effect_t_rp TFEffectStrength(const effect_t_rp *r){
	if (r){ 
		return *r; 
	}
	double effect =getValueFromSysParam(system_parameters.rates.tf.effectSize);
	return effect;}

rate_t_rp messengerProductionRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.mess.production);
}
rate_t_rp microProductionRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.micro.production);
}
rate_t_rp proteinProductionRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.tf.production);
}
rate_t_rp messengerDecayRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.mess.decay);
}
rate_t_rp microDecayRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.micro.decay);
}
rate_t_rp proteinDecayRate(const rate_t_rp *r){
	return (r) ? *r : getValueFromSysParam(system_parameters.rates.tf.decay);
}

rate_t_rp miRAssocConst(const rate_t_rp *assoc){
	if (assoc)
		return *assoc;
	rate_t_rp scale = getValueFromSysParam(system_parameters.rates.micro.frequency);
	rate_t_rp aff = getValueFromSysParam(system_parameters.rates.micro.affinity);
	return (2.0*scale*aff)/(1.0+aff);
}

rate_t_rp miRDissocConst(const rate_t_rp * dissoc){
	if (dissoc)
		return *dissoc;
	rate_t_rp scale = getValueFromSysParam(system_parameters.rates.micro.frequency);
	rate_t_rp aff = getValueFromSysParam(system_parameters.rates.micro.affinity);
	return (2.0*scale)/(1.0+aff);
}

rate_t_rp TFAssocConst(const rate_t_rp *assoc){
	if (assoc)
		return *assoc;
	rate_t_rp scale = getValueFromSysParam(system_parameters.rates.tf.frequency);
	rate_t_rp aff = getValueFromSysParam(system_parameters.rates.tf.affinity);
	return (2.0*scale*aff)/(1.0+aff);
}

rate_t_rp TFDissocConst(const rate_t_rp * dissoc){
	if (dissoc)
		return *dissoc;
	rate_t_rp scale = getValueFromSysParam(system_parameters.rates.tf.frequency);
	rate_t_rp aff = getValueFromSysParam(system_parameters.rates.tf.affinity);
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
    getRandomDiscrete(getStream(),nElems,parms->connectionDistribution,toRet.data);
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

void ConstOrCtsDistribution_Free(ConstOrCtsDistribution ** toFree){
    ConstOrCtsDistribution * loc = *toFree;
    if(loc){

    if(loc->isConst == -1){
        if(loc->distn){
            CompoundContinuousDistribution_free(loc->distn);
            free(loc->distn);
        }
        loc->distn = NULL;        
    }else if (loc->isConst == 1){
        loc->rate = 0.0;

    }else{
        fprintf(stderr,"trying to free compoundcts which is not initialized\n");
        exit(45436);
    }
        loc->isConst = 0;
        free(loc);
    }
    *toFree = NULL;
}


static void ConnectionParametersUnit_free(struct ConnectionParametersUnit **toFree){
    if(*toFree){
    if((*toFree)->connectionDistribution)
	{
        CompoundDiscreteDistribution_free((*toFree)->connectionDistribution);
        free((*toFree)->connectionDistribution);
        (*toFree)->connectionDistribution = NULL;
    }
	(*toFree)->connTp = UNDEFINED;
    free(*toFree);
    *toFree = NULL;
    }
}
static void connectionParameters_free(connectionParameters *toFree){
    ConnectionParametersUnit_free(&toFree->miRConns);
    ConnectionParametersUnit_free(&toFree->tfCodingConns);
    ConnectionParametersUnit_free(&toFree->tfNoncodingConns);
	
}

static void TFRates_free(struct TFRates * toFree){
     ConstOrCtsDistribution_Free(&toFree->production);
    ConstOrCtsDistribution_Free(&toFree->decay);
    ConstOrCtsDistribution_Free(&toFree->affinity);
    ConstOrCtsDistribution_Free(&toFree->frequency);
    ConstOrCtsDistribution_Free(&toFree->effectSize);

}

static void MicroRates_free(struct MicroRates * toFree){
     ConstOrCtsDistribution_Free(&toFree->production);
    ConstOrCtsDistribution_Free(&toFree->decay);
    ConstOrCtsDistribution_Free(&toFree->affinity);
    ConstOrCtsDistribution_Free(&toFree->frequency);
    ConstOrCtsDistribution_Free(&toFree->prodEffectStrength);
    ConstOrCtsDistribution_Free(&toFree->decayEffectStrength);
    //toFree->production =toFree->decay =toFree->affinity =toFree->frequency =
      //  toFree->prodEffectStrength =toFree->decayEffectStrength = NULL;
}

static void MessRates_free(struct MessRates * toFree){
    ConstOrCtsDistribution_Free(&toFree->production);
    ConstOrCtsDistribution_Free(&toFree->decay);
    //toFree->production = toFree->decay = NULL;
}

static void defaultRates_free(struct defaultRates * toFree){
    MessRates_free(&toFree->mess);
    MicroRates_free(&toFree->micro);
    TFRates_free(&toFree->tf);
    ConstOrCtsDistribution_Free(&toFree->defaultRate);
}


void destroyParameters() {
	connectionParameters_free(&system_parameters.connParms);
    defaultRates_free(&system_parameters.rates);
    memset(&system_parameters.initialQtys,0,sizeof(system_parameters.initialQtys));
}
