#ifndef NETWORK_DISTRIBUTIONS_RP__ 
#define NETWORK_DISTRIBUTIONS_RP__
#include <mkl_vsl.h>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include "../precision.h"
struct myContinuousDistributionArgs;
struct myDiscreteDistributionArgs;

typedef struct CompoundContinuousDistribution{
    struct myContinuousDistributionArgs * const distns;
    numeric_t_rp * const weights;
    const numeric_t_rp sumWeights;
    const int len;
} CompoundContinuousDistribution;

typedef struct CompoundDiscreteDistribution{
   struct  myDiscreteDistributionArgs * distns;
    numeric_t_rp * weights;
    const numeric_t_rp sumWeights;
    const int len;
} CompoundDiscreteDistribution;

//void parseXMLDiscreteDistribution(xmlDocPtr doc, xmlNodePtr cur,myDiscreteDistributionArgs * myDiscArgs);
//void parseXMLContinuousDistribution(xmlDocPtr doc, xmlNodePtr cur,myContinuousDistributionArgs *myCtsArgs);

CompoundContinuousDistribution *parseXMLCompoundContinuousDistribution(xmlDocPtr doc,xmlNodePtr cur,const char * parentName);
CompoundDiscreteDistribution *parseXMLCompoundDiscreteDistribution(xmlDocPtr doc,xmlNodePtr cur);
void getRandomDiscrete(VSLStreamStatePtr stream,int  numElts, CompoundDiscreteDistribution * distn, int * res);
void getRandomContinuous(VSLStreamStatePtr stream,int  numElts, CompoundContinuousDistribution * distn, numeric_t_rp * res);
void CompoundDiscreteDistribution_free(CompoundDiscreteDistribution * toFree);
void CompoundContinuousDistribution_free(CompoundContinuousDistribution * toFree);


#endif