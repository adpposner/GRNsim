#ifndef PARAMETERS_H_RP__
#define PARAMETERS_H_RP__
#include "globals.h"
#include <stdio.h>



typedef struct genesDims {
	size_t nDNA,nMess,nMicro,nProt;
} genesDims;







extern genesDims globalDims;

extern int networkGenSeed;


qty_type getInitialCodingQty();
qty_type getInitialNonCodingQty();
qty_type getInitialMessQty();
qty_type getInitialMicroQty();
qty_type getInitialQty(species_t sp);

effect_t_rp miRProdEffectStrength(const effect_t_rp *r); 
effect_t_rp miRDecayEffectStrength(const effect_t_rp *r); 
rate_t_rp TFAffinity(const rate_t_rp *r); 
effect_t_rp TFEffectStrength(const effect_t_rp *r); 


rate_t_rp messengerProductionRate(const rate_t_rp *r);
rate_t_rp microProductionRate(const rate_t_rp *r);
rate_t_rp proteinProductionRate(const rate_t_rp *r);
rate_t_rp messengerDecayRate(const rate_t_rp *r);
rate_t_rp microDecayRate(const rate_t_rp *r);
rate_t_rp proteinDecayRate(const rate_t_rp *r);

rate_t_rp miRAssocConst(const rate_t_rp *assoc);
rate_t_rp miRDissocConst(const rate_t_rp * dissoc);
rate_t_rp TFAssocConst(const rate_t_rp *assoc);
rate_t_rp TFDissocConst(const rate_t_rp * dissoc);


unsigned char connsByOutdegree(species_t sp);

UnsignedIntArray getRandomConnectionNumbersForElements(unsigned int leftLen, unsigned int rightLen, species_t left,species_t right);
void initializeParameters(const char * cwd);
void destroyParameters();
unsigned char * networkParamHash(unsigned char includeNetGenSeed);

#endif