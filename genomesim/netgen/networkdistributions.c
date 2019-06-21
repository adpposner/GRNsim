#include <string.h>
#include <mkl_vsl.h>
#include <libxml/tree.h>
#include <libxml/parser.h>
#include "../include/globals.h"
#include "../include/netgen/networkdistributions.h"
#include "../include/randomnumbers.h"
#include <mathimf.h>

enum DistributionTypeEnum{
    CONTINUOUS_DISTRIBUTION_CDE=0,DISCRETE_DISTRIBUTION_DDE=1
};

typedef enum ContinuousDistributionEnum{
    CONSTANT_CDE=0, UNIFORM_CDE, GAUSSIAN_CDE,EXPONENTIAL_CDE,
    LOGNORMAL_CDE,GAMMA_CDE,BETA_CDE,UNKNOWN_CDE
} ContinuousDistributionEnum;


typedef enum DiscreteDistributionEnum{
    CONSTANT_DDE=0, UNIFORM_DDE,BERNOULLI_DDE,
    BINOMIAL_DDE,POISSON_DDE,NEGATIVE_BINOMIAL_DDE,UNKNOWN_DDE
} DiscreteDistributionEnum;

const char * discrete_distribution_names[] = {"constant","uniform","bernoulli","binomial","poisson","negbinomial",NULL};
const char * continuous_distribution_names[] = {"constant","uniform","gaussian","exponential","lognormal","gamma","beta",NULL};



// int verifyDistributionName(xmlChar * distnname,DiscreteDistributionEnum expected){
//     int toRet = 1;
//     toRet = xmlStrcmp(distnname,(const xmlChar *)discrete_distribution_names[expected]);
//     if (!toRet)
//         return 0;
//     else{
//         fprintf(stderr, )
//     }
// };

#define MAX_DISCRETE_INT_ARGS 3
#define MAX_DISCRETE_DBL_ARGS 2
#define MAX_CTS_DBL_ARGS 5

const char * const constIntArgNames[3] ={"value",NULL,NULL};
const char * const  unifIntArgNames[3] ={"minVal","maxVal",NULL};
const char * const  bernoulliIntArgNames[3] ={"valFail","valSuccess",NULL};
const char * const binomialIntArgNames[3] ={"minVal","maxVal",NULL};
const char * const  poissonIntArgNames[3] ={"minVal","maxVal",NULL};
const char * const negBinomialIntArgNames[4] = {"minVal","maxVal","nSuccess",NULL};

const char * const constDblArgNames[2] ={NULL,NULL};
const char * const  unifDblArgNames[2] ={NULL,NULL};
const char * const  bernoulliDblArgNames[2] ={"successP",NULL};
const char * const binomialDblArgNames[2] ={"successP",NULL};
const char * const  poissonDblArgNames[2] ={"lambda",NULL};
const char * const negBinomialDblArgNames[2] = {"successP",NULL};
const char * const * discIntArgNames_[] = {constIntArgNames,unifIntArgNames,bernoulliIntArgNames,binomialIntArgNames,poissonIntArgNames,negBinomialIntArgNames};
const char * const * discDblArgNames_[] = {constDblArgNames,unifDblArgNames,bernoulliDblArgNames,binomialDblArgNames,poissonDblArgNames,negBinomialDblArgNames};

//const char * continuousIntArgNames_[][1] = {{NULL},{NULL},{NULL},{NULL},{NULL},{NULL},{NULL}};

    const char * const constCtsDblArgNames[5] = {"value",NULL,NULL,NULL,NULL};
const char * const unifCtsDblArgNames[5] = {"minVal","maxVal",NULL,NULL,NULL};

const char * const gaussianCtsDblArgNames[5] = {"mean","sd",NULL,NULL,NULL};
const char * const expCtsDblArgNames[5] = {"displacement","scale",NULL,NULL,NULL};
const char * const lgnormalCtsDblArgNames[5] = {"mu","sigma","displacement","scale",NULL};

const char * const gammaCtsDblArgNames[5] = {"shape","displacement","scale",NULL,NULL};
const char * const betaCtsDblArgNames[5] = {"shapea","shapeb","displacement","scale",NULL};

const char * const * continuousDblArgNames_[] = {constCtsDblArgNames,unifCtsDblArgNames,gaussianCtsDblArgNames,expCtsDblArgNames,
lgnormalCtsDblArgNames,gammaCtsDblArgNames,betaCtsDblArgNames};

const numeric_t_rp defaultConstCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp defaultUnifCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp defaultGaussianCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp defaultExpCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp defaultLgnormalCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,0.0,INFINITY,INFINITY};
const numeric_t_rp defaultGammaCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp defaultBetaCtsDblArgVals[MAX_CTS_DBL_ARGS] = {INFINITY,INFINITY,INFINITY,INFINITY,INFINITY};
const numeric_t_rp * const  defaultCtsDblArgVals_[] =  {defaultConstCtsDblArgVals, defaultUnifCtsDblArgVals, 
defaultGaussianCtsDblArgVals, defaultExpCtsDblArgVals, defaultLgnormalCtsDblArgVals, defaultGammaCtsDblArgVals,
defaultBetaCtsDblArgVals};

#define MAXINTORDBLARGS     10
typedef struct myDiscreteDistributionArgs{
    const DiscreteDistributionEnum myDistType;
    const char * const * intArgNames;
    const int intArgs[MAX_DISCRETE_INT_ARGS];
    const char * const * dblArgNames;
    const numeric_t_rp dblArgs[MAX_DISCRETE_DBL_ARGS];
} myDiscreteDistributionArgs;

typedef struct myContinuousDistributionArgs{
    const ContinuousDistributionEnum myDistType;
    const char * const * dblArgNames;
    const numeric_t_rp dblArgs[MAX_CTS_DBL_ARGS];
} myContinuousDistributionArgs;



int getIndexOfArgName(const char * argName, const char * const * argslist){
    int i;
    const char * const * currArg = argslist;
    for(i=0;i<MAXINTORDBLARGS;i++,currArg++){
        if(!(*currArg)){
            fprintf(stderr,"%s:%d - could not find argument %s in argslist\nValid arguments include:",__FILE__,__LINE__,argName);
            for(currArg = argslist;*currArg;currArg++){
                fprintf(stderr,"%s, ",*currArg);
            }
            fputc('\n',stderr);
            exit(-245);
        }else if(!strcmp(argName,*currArg))
                break;
    }
    return i;
}

static numeric_t_rp getCtsDefaultArgument(ContinuousDistributionEnum distType,const char * argName){
    int i;
    i=getIndexOfArgName(argName,continuousDblArgNames_[distType]);
    numeric_t_rp toRet;
    toRet = defaultCtsDblArgVals_[distType][i];
    if(isfinite(toRet))
        fprintf(stdout,"Argument %s for continuous distribution %s not specified - default value found = %f\n",
        argName,continuous_distribution_names[distType],toRet);
    return toRet;
}


int getDiscreteIntegerArgument(const char * argName,myDiscreteDistributionArgs * args){
    int i;
    i= getIndexOfArgName(argName,args->intArgNames);
    return args->intArgs[i];
}

numeric_t_rp getDiscreteNumericArgument(const char * argName,myDiscreteDistributionArgs * args){
    int i;
    i=getIndexOfArgName(argName,args->dblArgNames);
    return args->dblArgs[i];
}

int setDiscreteIntegerArgument(const char * argName,int value, myDiscreteDistributionArgs * args){
    int i;
    i= getIndexOfArgName(argName,args->intArgNames);
    *(int *)(&args->intArgs[i]) = value;
    return i;
}

int setDiscreteNumericArgument(const char * argName,numeric_t_rp value, myDiscreteDistributionArgs * args){
    int i;
    i= getIndexOfArgName(argName,args->dblArgNames);
    *(numeric_t_rp *)(&args->dblArgs[i]) = value;
    return i;
}

numeric_t_rp getContinuousNumericArgument(const char * argName,myContinuousDistributionArgs * args){
    int i;
    i=getIndexOfArgName(argName,args->dblArgNames);
    return args->dblArgs[i];
}


int setContinuousNumericArgument(const char * argName,numeric_t_rp value, myContinuousDistributionArgs * args){
    int i;
    i= getIndexOfArgName(argName,args->dblArgNames);
    *(numeric_t_rp *)(&args->dblArgs[i]) = value;
    return i;
}



myDiscreteDistributionArgs emptyDiscreteDistributionArgs(const char * distnName){
    int i;
    const char ** currDistName = discrete_distribution_names;
    if(distnName){
    for(i=0; i<UNKNOWN_DDE;i++,currDistName++){
        if(!(*currDistName)){
            fprintf(stderr,"%s:%d - Distribution named %s not found\n",__FILE__,__LINE__, distnName);
            exit(-3425);
        }
        else if(!strcmp(distnName,*currDistName))
            break;
    }
    myDiscreteDistributionArgs toRet = {.myDistType = i,.intArgs={0},.intArgNames=discIntArgNames_[i],.intArgs={0},.dblArgNames=discDblArgNames_[i],.dblArgs={0.0}};
    return toRet;}
    else
    {
        myDiscreteDistributionArgs toRet = {.myDistType = UNKNOWN_DDE,
        .intArgNames=discIntArgNames_[0],
        .intArgs={0},.dblArgNames=discDblArgNames_[0],.dblArgs={0.0}};
    }
    
}


myContinuousDistributionArgs emptyContinuousDistributionArgs(const char * distnName){
    int i;
    const char ** currDistName = continuous_distribution_names;
    if(distnName){
    for(i=0;i<UNKNOWN_CDE;i++,currDistName++){
        if(!(*currDistName)){
            fprintf(stderr,"%s:%d -Continuous Distribution named %s not found\n",__FILE__,__LINE__, distnName);
            exit(-3428);
        }else if (!strcmp(distnName,*currDistName))
            break;
    }
    myContinuousDistributionArgs toRet = {.myDistType = i,.dblArgNames=continuousDblArgNames_[i],.dblArgs={0.0}};
    return toRet;}
    else{
            myContinuousDistributionArgs toRet = {.myDistType = UNKNOWN_CDE,.dblArgNames=continuousDblArgNames_[0],.dblArgs={0.0}};
    }
}

void printDiscreteDistributionArgs(FILE * fh,myDiscreteDistributionArgs distnargs){
    int i;
    i=distnargs.myDistType;
    int numPossDists = sizeof(discrete_distribution_names) / sizeof(discrete_distribution_names[0]);
    fprintf(fh,"sizeof discrete didtnnames = %d\n",numPossDists);
    if (i > numPossDists){
        fprintf(stderr,"invalid distribution type of value %d",i);
        exit(-232415);
    }
    fprintf(fh,"distribution of type: %s \n",discrete_distribution_names[i]);
    const char * const* iArgName;
    const int * currIntArg;
    for(currIntArg=distnargs.intArgs,iArgName = distnargs.intArgNames;;currIntArg++,iArgName++){
        if (!(*iArgName)){
            break;
        }
        fprintf(fh,"\t INTEGER: %s =%d\n",*iArgName,*currIntArg);
    }
    const char *const* dArgName;
    const numeric_t_rp * currDblArg;
    for(currDblArg=distnargs.dblArgs,dArgName = distnargs.dblArgNames;;currDblArg++,dArgName++){
        if (!(*dArgName)){
            break;
        }
        fprintf(fh,"\t numeric_t_rp: %s =%f\n",*dArgName,*currDblArg);
    }
}

void printContinuousDistributionArgs(FILE * fh,myContinuousDistributionArgs distnargs){
    int i;
    i=distnargs.myDistType;
    int numPossDists = sizeof(continuous_distribution_names) / sizeof(continuous_distribution_names[0]);
    fprintf(fh,"sizeof cts didtnnames = %d\n",numPossDists);
    if (i > numPossDists){
        fprintf(stderr,"invalid distribution type of value %d",i);
        exit(-232415);
    }
    fprintf(fh,"distribution of type: %s \n",continuous_distribution_names[i]);
    const char *const * dArgName;
    const numeric_t_rp * currDblArg;
    for(currDblArg=distnargs.dblArgs,dArgName = distnargs.dblArgNames;;currDblArg++,dArgName++){
        if (!(*dArgName)){
            break;
        }
        fprintf(fh,"\t numeric_t_rp: %s =%f\n",*dArgName,*currDblArg);
    }
}


CompoundDiscreteDistribution * CompoundDiscreteDistribution_alloc(int len){
    
    CompoundDiscreteDistribution * toRet = calloc(1,sizeof(*toRet));


    assert(toRet);
    //printf("sizeof discargs=%zu,sizeof struct = %zu\n",sizeof(*discArgs),sizeof(myDiscreteDistributionArgs));

    DROPCONST(numeric_t_rp *,toRet->weights) = calloc(len,sizeof(*toRet->weights));
    DROPCONST(myDiscreteDistributionArgs *,toRet->distns) = calloc(len,sizeof(*toRet->distns));
    assert(toRet->distns);
    assert(toRet->weights);
    DROPCONST(numeric_t_rp,toRet->sumWeights)=0.0;
    DROPCONST(int,toRet->len) = len;
    myDiscreteDistributionArgs emp = emptyDiscreteDistributionArgs(NULL);
    for(--len;len>=0;len--)
         memcpy(toRet->distns+len,&emp,sizeof(emp));
    return toRet;
}

CompoundContinuousDistribution *CompoundContinuousDistribution_alloc(int len){
    myContinuousDistributionArgs * ctsArgs = calloc(len,sizeof(*ctsArgs));
    numeric_t_rp * myWeights = calloc(len,sizeof(*myWeights));
    CompoundContinuousDistribution *toRet = calloc(1,sizeof(*toRet));
    DROPCONST(numeric_t_rp *,toRet->weights) = myWeights;
    DROPCONST(myContinuousDistributionArgs *,toRet->distns) = ctsArgs;
    DROPCONST(int,toRet->len) = len;
    DROPCONST(numeric_t_rp,toRet->sumWeights) = 0.0;
    myContinuousDistributionArgs emp = emptyContinuousDistributionArgs(NULL);
    for(--len;len>=0;len--)
        memcpy(toRet->distns+len,&emp,sizeof(emp));
    return toRet;
}

void CompoundDiscreteDistribution_free(CompoundDiscreteDistribution * toFree){
    if(toFree->weights)
        free(toFree->weights);
    if (toFree->distns)
        free(toFree->distns);
    DROPCONST(numeric_t_rp *,toFree->weights) = NULL;
    DROPCONST(myContinuousDistributionArgs *, toFree->distns) = NULL;
    DROPCONST(int, toFree->len) = 0;
    DROPCONST(numeric_t_rp, toFree->sumWeights) = 0.0;
}

void CompoundContinuousDistribution_free(CompoundContinuousDistribution * toFree){
    if(toFree->weights)
        free(toFree->weights);
    if (toFree->distns)
        free(toFree->distns);
    DROPCONST(numeric_t_rp *,toFree->weights) = NULL;
    DROPCONST(myContinuousDistributionArgs *,toFree->distns) = NULL;
    DROPCONST(int, toFree->len) = 0;
    DROPCONST(numeric_t_rp, toFree->sumWeights) = 0.0;
}

static int getXMLintAttribute(xmlDocPtr doc, xmlNodePtr cur, const char * attribName){
    xmlChar *key;
    int toRet=0;
    key = xmlGetProp(cur,attribName);
    if(key){
        toRet = strtol((const char *) key,NULL,10);
        xmlFree(key);
    }else{
        return -1;
    }
    return toRet;
}

static int XMLlistAttributes(xmlNodePtr cur){
    xmlAttr * attri = cur->properties;
    while(attri){
        if (attri->name && attri->children){
            xmlChar * attrVal;
            attrVal = xmlNodeListGetString(cur->doc,attri->children,1);
            if(attrVal)
                {printf("attribute %s = %s\n",attri->name,attrVal);xmlFree(attrVal);}
        }
        attri=attri->next;
    }
}

static numeric_t_rp getXMLNumericAttribute(xmlDocPtr doc, xmlNodePtr cur, const char * attribName){
    xmlChar *key;
    numeric_t_rp toRet=0;
    key = xmlGetProp(cur,attribName);
    if(key){

			toRet = strtod((const char *)key, NULL);

			xmlFree(key);
    }else{
        return -1;
    }
    
    return toRet;
}


static int isValidDistributionName(enum DistributionTypeEnum distnType,const char * dname){
    const char ** origDistn = (distnType == CONTINUOUS_DISTRIBUTION_CDE) ? continuous_distribution_names : discrete_distribution_names;
    const char ** validDistns = origDistn;
    int i;
    for(i=0;*validDistns;validDistns++,i++){
        if(!strcmp(*validDistns,dname)){
            break;
        }
    }
    if(*validDistns)
        return i;
    else
    {
        fprintf(stderr,"%s name for %s distribution type invalid. Available types are: ",dname,(distnType == CONTINUOUS_DISTRIBUTION_CDE)? "continuous":"discrete");
        for(validDistns=origDistn;*validDistns;validDistns++){
            fprintf(stderr, "%s, ",*validDistns);
        }
        fputc('\n',stderr);
        exit(-62839);
    }
    
}


void printDiscreteDistArgNames(FILE * fh, const char * distname){
    int i,j;
    i=isValidDistributionName(DISCRETE_DISTRIBUTION_DDE,distname);
    fprintf(fh,"Integer arguments for %s distribution: ",distname);
    const char *const* argslist = discIntArgNames_[i];
    for(;*argslist;argslist++){
        fprintf(fh,"%s, ",*argslist);
    }
    fputc('\n',fh);
    fprintf(fh,"numeric_t_rp arguments for %s distribution: ", distname);
    argslist = discDblArgNames_[i];
    for(;*argslist;argslist++){
        fprintf(fh,"%s, ",*argslist);
    }
    fputc('\n',fh);
}

void printContinuousDistArgNames(FILE * fh, const char * distname){
    int i,j;
    i=isValidDistributionName(CONTINUOUS_DISTRIBUTION_CDE,distname);
    fprintf(fh,"numeric_t_rp arguments for %s distribution: ", distname);
    const char * const * argslist;
    argslist = continuousDblArgNames_[i];
    for(;*argslist;argslist++){
        fprintf(fh,"%s, ",*argslist);
    }
    fputc('\n',fh);
}

static void getWeight(xmlDocPtr doc, xmlNodePtr cur, numeric_t_rp * weightVal){
    numeric_t_rp wv;
    wv = getXMLNumericAttribute(doc,cur,"weight");
    if (wv < 0.0){
    
        fprintf(stderr,"must specify weights for all disctributions,%s %f\n",(const char *) cur->name,wv);
        exit(-12458);
    }
    *weightVal = wv;
}

void parseXMLDiscreteDistribution(xmlDocPtr doc, const xmlNodePtr cur,myDiscreteDistributionArgs * myDiscArgs){
    
    //XMLlistAttributes(cur);
    assert(!xmlStrcmp(cur->name,(const xmlChar *)"distribution"));
    xmlChar * key;
    key = xmlGetProp(cur,"type");
    myDiscreteDistributionArgs dArgs = emptyDiscreteDistributionArgs(NULL);
         

    if(!key){
        fprintf(stderr,"type attribute not found in distribution XML\n");
        exit(-21322);
    }else{
        isValidDistributionName(DISCRETE_DISTRIBUTION_DDE,(const char *)key);
        myDiscreteDistributionArgs emp = emptyDiscreteDistributionArgs((const char *)key);
        memcpy(&dArgs,&emp,sizeof(dArgs));
        xmlFree(key);
    }
       

    const char * const * intArgName = dArgs.intArgNames;
    int iArgVal;
    for(;*intArgName;intArgName++){
        if ((iArgVal = getXMLintAttribute(doc,cur,*intArgName))<0){
            fprintf(stderr,"Integer argument %s not found for discrete distribution in xml attributes",*intArgName);
            printDiscreteDistArgNames(stderr,discrete_distribution_names[dArgs.myDistType]);
            exit(-48245);
        }else{
            setDiscreteIntegerArgument(*intArgName,iArgVal,&dArgs);
        }
    }
       
    numeric_t_rp dArgVal;
    const char * const * dblArgName = dArgs.dblArgNames;
    for(;*dblArgName;dblArgName++){
        if((dArgVal = getXMLNumericAttribute(doc,cur,*dblArgName))<0.0){
            fprintf(stderr,"numeric_t_rp argument %s not found for discrete distribution in xml attributes",*dblArgName);
            printDiscreteDistArgNames(stderr,discrete_distribution_names[dArgs.myDistType]);
            exit(-423438);
        }else{
            setDiscreteNumericArgument(*dblArgName,dArgVal,&dArgs);
        }
    }
       
        memmove(myDiscArgs,&dArgs,sizeof(dArgs));
       
}

void parseXMLContinuousDistribution(xmlDocPtr doc, xmlNodePtr cur,myContinuousDistributionArgs *myCtsArgs){
    assert(!xmlStrcmp(cur->name,(const xmlChar *)"distribution"));
    xmlChar * key;
    key = xmlGetProp(cur,"type");
    myContinuousDistributionArgs cArgs = emptyContinuousDistributionArgs(NULL);
    if(!key){
        fprintf(stderr,"type attribute not found in distribution XML\n");
        exit(-21322);
    }else{
        isValidDistributionName(CONTINUOUS_DISTRIBUTION_CDE,(const char *)key);
        myContinuousDistributionArgs emp = emptyContinuousDistributionArgs((const char *)key);
        memcpy(&cArgs,&emp,sizeof(cArgs));
        xmlFree(key);
    }
    numeric_t_rp dArgVal;
    const char * const * dblArgName = cArgs.dblArgNames;

    for(;*dblArgName;dblArgName++){
        if((dArgVal = getXMLNumericAttribute(doc,cur,*dblArgName))<0.0){
            if(isfinite(dArgVal = getCtsDefaultArgument(cArgs.myDistType,*dblArgName)))
                setContinuousNumericArgument(*dblArgName,dArgVal,&cArgs);
            else{
            fprintf(stderr,"numeric_t_rp argument %s not found for discrete distribution in xml attributes",*dblArgName);
            printContinuousDistArgNames(stderr,continuous_distribution_names[cArgs.myDistType]);
            exit(-4234238);
            }
        }else{
            setContinuousNumericArgument(*dblArgName,dArgVal,&cArgs);
        }
    }
    memcpy(myCtsArgs,&cArgs,sizeof(cArgs));
}

void setWeightSumForCompoundContinuousDistribution(CompoundContinuousDistribution *toSet){
    numeric_t_rp weightsum=0;
    int i;
    for(i=0;i<toSet->len;i++){
        weightsum+= toSet->weights[i];
    }
    *(numeric_t_rp *)&toSet->sumWeights = weightsum;
}

void setWeightSumForCompoundDiscreteDistribution(CompoundDiscreteDistribution *toSet){
    numeric_t_rp weightsum=0;
    int i;
    for(i=0;i<toSet->len;i++){
        weightsum+= toSet->weights[i];
    }
    *(numeric_t_rp *)&toSet->sumWeights = weightsum;
}


CompoundDiscreteDistribution * parseXMLCompoundDiscreteDistribution(xmlDocPtr doc,xmlNodePtr cur){
    while(cur!=NULL){
        if(!xmlStrcmp(cur->name,(const xmlChar * )"discrete_distribution"))
            break;
        cur=cur->next;
    }
    if(!cur){
        fprintf(stderr,"No discrete distribution found in config file!\n");
        exit(-28283);
    }
    int numElts = xmlChildElementCount(cur);
    //fprintf(stdout,"discrete compound distribution found, numelts=%d\n",numElts);
    CompoundDiscreteDistribution *toRet = CompoundDiscreteDistribution_alloc(numElts);
    xmlNodePtr curChild;
    int i;
    myDiscreteDistributionArgs * curDistn = toRet->distns;
    numeric_t_rp * curWeight = toRet->weights;
    *curWeight = 2.0;
    for(curChild=cur->children,i=0;curChild;curChild=curChild->next){
        if(!xmlStrcmp(curChild->name,(const xmlChar *)"distribution")){
            //printf("%s:%d - nodename = %s\n",__FILE__,__LINE__,(const char *)curChild->name);
            getWeight(doc,curChild,curWeight);
            //printf("%s:%d - nodename = %s\n",__FILE__,__LINE__,(const char *)curChild->name);
            const xmlNodePtr cc = curChild;
            parseXMLDiscreteDistribution(doc,cc,curDistn);
            assert(cc == curChild);
            //printf("%s:%d - nodename = %s,weight=%f\n",__FILE__,__LINE__,curChild->name,*curWeight);
            //printDiscreteDistributionArgs(stdout,*curDistn);
            curDistn++,curWeight++;
            i++;
            if (i==numElts) break;
        }
    }
    if(i != numElts){
        fprintf(stderr,"Invalid length in discrete distribution parsing numElts=%d i = %d \n",numElts,i);
        exit(-23352);
    }
    setWeightSumForCompoundDiscreteDistribution(toRet);
    return toRet;
}

CompoundContinuousDistribution * parseXMLCompoundContinuousDistribution(xmlDocPtr doc,xmlNodePtr cur,const char * parentName){
    while(cur!=NULL){
        if(!xmlStrcmp(cur->name,(const xmlChar * )"continuous_distribution"))
            break;
        cur=cur->next;
    }
    if(!cur){
        fprintf(stderr,"No continuous distribution found in config file for %s !\n",parentName);
        exit(-282383);
    }
    int numElts = xmlChildElementCount(cur);
    
    CompoundContinuousDistribution * toRet =
    CompoundContinuousDistribution_alloc(numElts);
    xmlNodePtr curChild;
    int i;
    myContinuousDistributionArgs * curDistn = toRet->distns;
    numeric_t_rp * curWeight = toRet->weights;
    for(curChild=cur->children,i=0;curChild;curChild=curChild->next){
        if(!xmlStrcmp(curChild->name,(const xmlChar *)"distribution")){
            //printf("%s:%d - nodename = %s\n",__FILE__,__LINE__,(const char *)curChild->name);
            getWeight(doc,curChild,curWeight);
            //printf("%s:%d - nodename = %s\n",__FILE__,__LINE__,(const char *)curChild->name);
            const xmlNodePtr cc = curChild;
            parseXMLContinuousDistribution(doc,cc,curDistn);
            assert(cc == curChild);
            //printf("%s:%d - nodename = %s,weight=%f\n",__FILE__,__LINE__,curChild->name,*curWeight);
            //printDiscreteDistributionArgs(stdout,*curDistn);
            curDistn++,curWeight++;
            i++;
            if (i==numElts) break;
        }
    }
    if(i != numElts){
        fprintf(stderr,"Invalid length in continuous distribution parsing\n");
        exit(-23352);
    }
    setWeightSumForCompoundContinuousDistribution(toRet);
    return toRet;
}






// int main(void){
//     myDiscreteDistributionArgs p = emptyDiscreteDistributionArgs("poisson");
//     setDiscreteNumericArgument("lambda",2.5,&p);
//     printDiscreteDistributionArgs(stdout,p);
//     myContinuousDistributionArgs cp = emptyContinuousDistributionArgs("beta");
//     setContinuousNumericArgument("scalefactor",2.5,&cp);
//     printContinuousDistributionArgs(stdout,cp);
//     return 0;
// }




int callConstantDiscrete(VSLStreamStatePtr stream, myDiscreteDistributionArgs * data, int * res){
    *res = getDiscreteIntegerArgument("value",data);
    return *res;
}

int callUniformDiscrete(VSLStreamStatePtr stream, myDiscreteDistributionArgs *data, int* res){
    assert(data->myDistType == UNIFORM_DDE);
    return viRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,res,
    getDiscreteIntegerArgument("minVal",data),getDiscreteIntegerArgument("maxVal",data));
}

int callBernoulliDiscrete(VSLStreamStatePtr stream,myDiscreteDistributionArgs * data,int * res){
    int status;
    status = viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF,stream,1,res,getDiscreteNumericArgument("successP",data));
    int successVal,failVal;
    successVal = getDiscreteIntegerArgument("valSuccess",data);
    failVal = getDiscreteIntegerArgument("valFail",data);
    if(*res){
        *res = (successVal >= 0) ? successVal : 1;
    }else{
        *res = (failVal >= 0) ? failVal : 0;
    }
    return status;
}


int callNegBinomialDiscrete(VSLStreamStatePtr stream, myDiscreteDistributionArgs * data, int * res){
    int status,tmp,minVal,maxVal;
    status = viRngNegBinomial(VSL_RNG_METHOD_NEGBINOMIAL_NBAR,stream,1,&tmp,getDiscreteIntegerArgument("nSuccess",data),
    getDiscreteNumericArgument("successP",data));
    minVal = getDiscreteIntegerArgument("minVal",data);
    maxVal = getDiscreteIntegerArgument("maxVal",data);
    tmp += minVal;
    tmp = (tmp > maxVal) ? maxVal : tmp;
    *res = tmp;
    return status;
}

int callBinomialDiscrete(VSLStreamStatePtr stream, myDiscreteDistributionArgs * data, int * res){
    int status;
    int tmp;
    int valDelta;
    valDelta = getDiscreteIntegerArgument("maxVal",data) - getDiscreteIntegerArgument("minVal",data);
    status = viRngBinomial(VSL_RNG_METHOD_BINOMIAL_BTPE,stream,1,&tmp,valDelta,getDiscreteNumericArgument("successP",data));
    tmp += getDiscreteIntegerArgument("minVal",data);
    *res = tmp;
    return status;
}

int callPoissonDiscrete(VSLStreamStatePtr stream,myDiscreteDistributionArgs * data, int * res){
    int status;
    int tmp;
    status=viRngPoisson(VSL_RNG_METHOD_POISSON_PTPE,stream,1,&tmp,getDiscreteNumericArgument("lambda",data));
    tmp += getDiscreteIntegerArgument("minVal",data);
    int maxval = getDiscreteIntegerArgument("maxVal",data);
    tmp = (tmp > maxval) ? maxval : tmp;
    *res = tmp;
    return status;
}

static void getRandomIntFromDiscrete(VSLStreamStatePtr stream,myDiscreteDistributionArgs * data, int *res){
    switch(data->myDistType){
        case CONSTANT_DDE: callConstantDiscrete(stream,data,res);break;
        case UNIFORM_DDE: callUniformDiscrete(stream,data,res);break;
        case BERNOULLI_DDE: callBernoulliDiscrete(stream,data,res);break;
        case BINOMIAL_DDE: callBinomialDiscrete(stream,data,res);break;
        case POISSON_DDE: callPoissonDiscrete(stream,data,res);break;
        case NEGATIVE_BINOMIAL_DDE: callNegBinomialDiscrete(stream,data,res);break;
        default: fprintf(stderr,"%s:%d - Invalid discrete distribution value = %d\n",__FILE__,__LINE__,data->myDistType);
            exit(-754);
    }
}

static void getRandomDiscrete_element(VSLStreamStatePtr stream, CompoundDiscreteDistribution * compoundDist, int * res){
    numeric_t_rp sel;
    int status;
    int i;
    myDiscreteDistributionArgs *selectedDistribution = compoundDist->distns;
    status = vNumberRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,&sel,0.0,compoundDist->sumWeights);
    for(i=0;i<compoundDist->len;i++){
        if(sel > compoundDist->weights[i]){
            sel = sel - compoundDist->weights[i];
            selectedDistribution++;
        }else{
            break;
        }
    }
    if(i == compoundDist->len){
        fprintf(stderr,"Could not get a value to select from in distribution %f, sumweights %f\n",sel,compoundDist->sumWeights);
        exit(-22482);
    }
    getRandomIntFromDiscrete(stream,selectedDistribution,res);
}

void getRandomDiscrete(VSLStreamStatePtr stream,int  numElts, CompoundDiscreteDistribution * distn, int * res){
    for(;numElts>0;numElts--,res++){
        getRandomDiscrete_element(stream,distn,res);
    }
}


int callConstantContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    *res =  getContinuousNumericArgument("value",data);   
    return 0;
}

int callUniformContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    int status;
    status = vNumberRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,res,getContinuousNumericArgument("minVal",data),getContinuousNumericArgument("maxVal",data));
    return status;
}

int callGaussianContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    int status;
    status = vNumberRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF,stream,1,res,getContinuousNumericArgument("mean",data),getContinuousNumericArgument("sd",data));
    return status;
}

int callExponentialContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    int status;
    status = vNumberRngExponential(VSL_RNG_METHOD_EXPONENTIAL_ICDF,stream,1,res,getContinuousNumericArgument("displacement",data),getContinuousNumericArgument("scale",data));
    return status;
}

int callLgNormalContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    int status;
    status = vNumberRngLognormal(VSL_RNG_METHOD_LOGNORMAL_ICDF,stream,1,res,getContinuousNumericArgument("mu",data),
    getContinuousNumericArgument("sigma",data),
    getContinuousNumericArgument("displacement",data),
    getContinuousNumericArgument("scale",data));
    return status;
}
int callGammaContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
    int status;
    status = vNumberRngGamma(VSL_RNG_METHOD_GAMMA_GNORM,stream,1,res,getContinuousNumericArgument("shape",data),getContinuousNumericArgument("displacement",data),
    getContinuousNumericArgument("scale",data));
    return status;
}
int callBetaContinuous(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp * res){
        int status;
         status = vNumberRngBeta(VSL_RNG_METHOD_BETA_CJA,stream,1,res,getContinuousNumericArgument("shapea",data),getContinuousNumericArgument("shapeb",data),
                  getContinuousNumericArgument("displacement",data),getContinuousNumericArgument("scale",data));
        return status;
}

static void getRandomDoubleFromCts(VSLStreamStatePtr stream,myContinuousDistributionArgs * data, numeric_t_rp *res){
    switch(data->myDistType){
        case CONSTANT_CDE: callConstantContinuous(stream,data,res);break;
        case UNIFORM_CDE: callUniformContinuous(stream,data,res);break;
        case GAUSSIAN_CDE: callGaussianContinuous(stream,data,res);break;
        case EXPONENTIAL_CDE: callExponentialContinuous(stream,data,res);break;
        case LOGNORMAL_CDE: callLgNormalContinuous(stream,data,res);break;
        case GAMMA_CDE: callGammaContinuous(stream,data,res);break;
        case BETA_CDE: callBetaContinuous(stream,data,res);break;
        default: fprintf(stderr,"%s:%d - Invalid continuous distribution value = %d\n",__FILE__,__LINE__,data->myDistType);
            exit(-7545);
    }
}

static void getRandomContinuous_element(VSLStreamStatePtr stream, CompoundContinuousDistribution * compoundDist, numeric_t_rp *res){
    int status;
    int i;
    numeric_t_rp sel;
    myContinuousDistributionArgs *selectedDistribution = compoundDist->distns;
    status = vNumberRngUniform(VSL_RNG_METHOD_UNIFORM_STD,stream,1,&sel,0.0,compoundDist->sumWeights);
    for(i=0;i<compoundDist->len;i++){
        if(sel > compoundDist->weights[i]){
            sel = sel - compoundDist->weights[i];
            selectedDistribution++;
        }else{
            break;
        }
    }
    if(i == compoundDist->len){
        fprintf(stderr,"Could not get a value to select from in distribution %f, sumweights %f\n",sel,compoundDist->sumWeights);
        exit(-2282);
    }
    getRandomDoubleFromCts(stream,selectedDistribution,res);
}

void getRandomContinuous(VSLStreamStatePtr stream,int  numElts, CompoundContinuousDistribution * distn, numeric_t_rp * res){
    for(;numElts>0;numElts--,res++){
        getRandomContinuous_element(stream,distn,res);
    }
}



// int callGeometricDiscrete(VSLStreamStatePtr stream,DiscreteDistribution_rp data, int * res){
//     int status;
//     int tmp;
//     status = viRngGeometric(VSL_RNG_METHOD_GEOMETRIC_ICDF,stream,1,&tmp,data.geometricArgs.successP);
//     if (data.geometricArgs.minVal > 0)
//         tmp = tmp + data.geometricArgs.minVal;
//     if (data.geometricArgs.maxVal > 0)
//         tmp = (tmp > data.geometricArgs.maxVal) ? data.geometricArgs.maxVal : tmp;
//     *res = tmp;
//     return status;
// }



// int callHyperGeometric(VSLStreamStatePtr stream,DiscreteDistribution_rp data, int * res){
//     int status;
//     viRngHypergeometric(VSL_RNG_METHOD_HYPERGEOMETRIC_H2PE,stream,1,res,data.hypergeomArgs.NtotalElements,data.hypergeomArgs.sampleSize,data.hypergeomArgs.numMarked);
//     return status;
// }




// typedef int (*DiscreteRNGFunc)(VSLStreamStatePtr,int,int,numeric_t_rp,numeric_t_rp,int,int*);

// static DiscreteRNGFunc rngFuncs[7] = {callUniformDiscrete,callBernoulliDiscrete,callGeometricDiscrete,callBinomial,callHyperGeometric,callPoissonDiscrete,callNegBinomialDiscrete};

// typedef struct CompositeDiscreteDistn_rp{
//     const int len;
//     DiscreteDistribution_rp * myDistributions;
//     numeric_t_rp * myWeights;
// } CompositeDiscreteDistn_rp;




// DiscreteDistribution_rp parseUniformDiscreteDistribution(xmlDocPtr doc, xmlNodePtr cur){
//     int minVal,maxVal;
//     if(verifyDistributionName(cur->name,UNIFORM_DDE)){
//         minVal = getXMLintAttribute(doc,cur,"minVal");
//         maxVal = getXMLintAttribute(doc,cur,"maxVal");
//     }
//     return makeUniformDiscreteDistribution(minVal,maxVal);
// }

// DiscreteDistribution_rp parseBernoulliDiscreteDistribution(xmlDocPtr doc, xmlNodePtr cur){
//     int valFail,valSuccess;
//     numeric_t_rp successP;
//     if(verifyDistributionName(cur->name,BERNOULLI_DDE)){
//         valFail = getXMLintAttribute(doc,cur,"valFail");
//         valSuccess = getXMLintAttribute(doc,cur,"valSuccess");
//         successP = getXMLNumericAttribute(doc,cur,"successP");
//     }
//     return makeBernoulliDistribution(valFail,valSuccess,successP);
// }

// DiscreteDistribution_rp parseGeometricDistribution(xmlDocPtr doc, xmlNodePtr cur){
//     int minVal,maxVal;
//     numeric_t_rp successP;
//     if(verifyDistributionName(cur->name,GEOMETRIC_DDE)){
//         minVal = getXMLintAttribute(doc,cur,"minVal");
//         maxVal = getXMLintAttribute(doc,cur,"maxVal");
//         successP = getXMLNumericAttribute(doc,cur,"successP");
//     }
//     return makeGeometricDistribution(minVal,maxVal,successP);
// }

// DiscreteDistribution_rp parseBinomialDistribution(xmlDocPtr doc, xmlNodePtr cur){
//     int minVal,maxVal;
//     numeric_t_rp successP;
//     if(verifyDistributionName(cur->name,BINOMIAL_DDE)){
//         minVal = getXMLintAttribute(doc,cur,"minVal");
//         maxVal = getXMLintAttribute(doc,cur,"maxVal");
//         successP = getXMLNumericAttribute(doc,cur,"successP");
//     }
//     return makeBinomialDistribution(minVal,maxVal,successP);
// }

// DiscreteDistribution_rp parsePoissonDistribution(xmlDocPtr doc, xmlNodePtr cur){
//     int minVal,maxVal;
//     numeric_t_rp lambda;
//     if(verifyDistributionName(cur->name,POISSON_DDE)){
//         minVal = getXMLintAttribute(doc,cur,"minVal");
//         maxVal = getXMLintAttribute(doc,cur,"maxVal");
//         lambda = getXMLNumericAttribute(doc,cur,"lambda");
//     }
//     return makePoissonDistribution(minVal,maxVal,lambda);
// }



// DiscreteDistribution_rp discreteDistributionFromXML(xmlDocPtr doc, xmlNodePtr cur){
//     if(!xmlStrcmp(cur->name,(const xmlChar *) "discrete_distribution")){

//     }else{
//         fprintf(stderr,"Failure in parsing discrete distribution, invalid class %s",(const char *) cur->name);
//     }
// }
