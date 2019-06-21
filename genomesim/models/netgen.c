#include "models/netgen.h"
#include "../include/models/commongeneticelement.h"
#include "parameters_p.h"
#include "randomnumbers.h"
#include "globals.h"
#include <string.h>

// #define VOIDCOMMONGENETICELEMENT  (CommonGeneticElement) { .name = {0}, .species = SPECIES_NONE, .qty = 0,          \
// .id = -1, .isEnabled = 0, .left = {.id = -1}, .right = {.id = -1}, .assocConstant = 0.0f, .dissocConstant = 0.0f,   \
//  .effectStrength = 0.0f, .decayConstant = 0.0f, .associatedIDs = NULL, .IDsAsPtrs = 0   };


typedef struct WrappedCGEA {
    CommonGeneticElementArray * myData;
    CommonGeneticElement * const codingStart;
    CommonGeneticElement * const nonCodingStart;
    CommonGeneticElement * const messStart;
    CommonGeneticElement * const mirStart;
    CommonGeneticElement * const tfStart;
    CommonGeneticElement * const end;
} WrappedCGEA;

CommonGeneticElementArray getArrayForSpecies(WrappedCGEA arr, species_t sp){
    switch(sp){
        case SPECIES_CODING: return (CommonGeneticElementArray){.data=arr.codingStart,.length = 
        arr.nonCodingStart - arr.codingStart};
        case SPECIES_NONCODING: return (CommonGeneticElementArray){.data=arr.nonCodingStart,.length = 
        arr.messStart - arr.nonCodingStart};
        case SPECIES_MESSENGER: return (CommonGeneticElementArray){.data=arr.messStart,.length = 
        arr.mirStart - arr.messStart};
        case SPECIES_MICRO: return (CommonGeneticElementArray){.data=arr.mirStart,.length = 
        arr.tfStart - arr.mirStart};
        case SPECIES_PROTEIN: return (CommonGeneticElementArray){.data = arr.tfStart,
        .length = arr.end - arr.tfStart};
        default: fprintf(stderr,"invalid %s:%d\n",__FILE__,__LINE__);
        assert(0);
    }
    return (CommonGeneticElementArray) {.data = NULL,.length = 0};
}

WrappedCGEA CGEA_calloc(const int nMess, const int nMicro){
    const int nTot = 3*nMess + 2*nMicro;
    CommonGeneticElementArray * tmparr = my_calloc(1,sizeof(*tmparr));
    CommonGeneticElementArray_calloc(tmparr,nTot);
    WrappedCGEA toRet = {.myData = tmparr,.codingStart= tmparr->data,
    .nonCodingStart = tmparr->data + nMess, .messStart = tmparr->data + nMess + nMicro,
    .mirStart = tmparr->data + 2*nMess + nMicro, .tfStart = tmparr->data + 2*nMess + 2*nMicro,
    .end = tmparr->data + 3*nMess + 2*nMicro};
    //set species, qty and names
    CommonGeneticElement * curr = res->data;
    int myID = 0;
    for(int k=0;k<nMess;k++,curr++){
        curr->id = myID++;
        curr->species = SPECIES_CODING;
        sprintf(curr->name,"gene %d",k);
        curr->produces.ptr = curr+nMess + nMicro;
        curr->IDsAsPtrs = 1;
    }
    for(int k=0;k<nMicro;k++,curr++){
        curr->id = myID++;
        curr->species = SPECIES_NONCODING;
        sprintf(curr->name,"ncgene %d",k);
        curr->produces.ptr = curr+nMess + nMicro;
        curr->IDsAsPtrs = 1;
    }
    for(int k=0;k<nMess;k++,curr++){
        curr->id = myID++;
        curr->species = SPECIES_MESSENGER;
        sprintf(curr->name,"mess %d",k);
        curr->produces.ptr = curr+nMess + nMicro;
        curr->producedBy.ptr = curr - nMess - nMicro;
        curr->IDsAsPtrs = 1;
    }
    for(int k=0;k<nMicro;k++,curr++){
        curr->id = myID++;
        curr->species = SPECIES_MICRO;
        sprintf(curr->name,"mir %d",k);
        curr->produces.ptr = NULL;
        curr->producedBy.ptr = curr - nMess - nMicro;
        curr->IDsAsPtrs = 1;
        
    }
    for(int k=0;k<nMess;k++,curr++){
        curr->id = myID++;
        curr->species = SPECIES_PROTEIN;
        sprintf(curr->name,"prot %d",k);
        curr->produces.ptr = NULL;
        curr->producedBy.ptr = curr - nMess - nMicro;
        curr->IDsAsPtrs = 1;
    }
    assert(curr == (res->data + res->length));
    return res;
}

static void fillRandomBinomialArray(struct ConnectionParametersUnit * parms, 
CommonGeneticElementArray * toFill){
    UnsignedIntArray *arr = malloc(sizeof(*arr));
    UnsignedIntArray_calloc(arr,toFill->length);
    randomBinomialArray(arr->length,parms->deg_low,parms->deg_high,parms->effect_prob,arr->data);
    unsigned * currRand = arr->data;
    ARRAYTY
        start->associatedIDs = my_calloc(1,sizeof*start->associatedIDs);
        CGE_IDPtrArray_calloc(start->associatedIDs,*currRand);
    UnsignedIntArray_free(arr);
    my_free(arr);
}

void setCGE_IDPtrArray_rdm(CGE_IDPtrArray * toSet, CGE_IDPtrArray * ptrsToShuff){
    shufflePtrArray(ptrsToShuff->data,ptrsToShuff->length,sizeof(*ptrsToShuff->data));
    memmove(toSet->data,ptrsToShuff->data,toSet->length * sizeof(*toSet->data));
}

CGE_IDPtrArray getPtrArrayFor(CommonGeneticElement * start, int length){
    CGE_IDPtrArray toRet;
    CGE_IDPtrArray_calloc(&toRet,length);
    for (int k=0;k<length;k++)
        toRet.data[k].ptr = start+k;
    return toRet;
}

void associatedIDs_calloc_rdm(CommonGeneticElementArray *toAlloc,const int nMess, const int nMicro){
 
    fillRandomBinomialArray(system_parameters.connParms.miRConns,toAlloc->data+nMess+nMicro,nMess);
    fillRandomBinomialArray(system_parameters.connParms.tfConns, toAlloc->data,nMess + nMicro);
    CGE_IDPtrArray myTFs =  getPtrArrayFor(toAlloc->data+2*nMess+2*nMicro,nMess);
    CGE_IDPtrArray myMirs = getPtrArrayFor(toAlloc->data + 2*nMess + nMicro, nMicro);
    for(int k=0;k<nMess+nMicro;k++){
        setCGE_IDPtrArray_rdm(toAlloc->data[k].associatedIDs,&myTFs);
    }
    for(int k=0;k<nMess;k++){
        setCGE_IDPtrArray_rdm(toAlloc->data[k].associatedIDs,&myMirs);
    }
    CGE_IDPtrArray_free(&myTFs);
    CGE_IDPtrArray_free(&myMirs);
}


