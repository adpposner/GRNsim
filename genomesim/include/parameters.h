//parameters.h - public interface for system parameters
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
#ifndef PARAMETERS_H_RP__
#define PARAMETERS_H_RP__
#include "globals.h"
#include <stdio.h>



typedef struct genesDims {
	size_t nDNA,nMess,nMicro,nProt;
} genesDims;







extern genesDims globalDims;

extern int networkGenSeed;


ulong_type getInitialCodingQty();
ulong_type getInitialNonCodingQty();
ulong_type getInitialMessQty();
ulong_type getInitialMicroQty();
ulong_type getInitialQty(species_t sp);

effect_t_rp miREffectStrength(effect_t_rp *r); 
rate_t_rp TFAffinity(rate_t_rp *r); 
effect_t_rp TFEffectStrength(effect_t_rp *r); 


rate_t_rp messengerProductionRate(rate_t_rp *r);
rate_t_rp microProductionRate(rate_t_rp *r);
rate_t_rp proteinProductionRate(rate_t_rp *r);
rate_t_rp messengerDecayRate(rate_t_rp *r);
rate_t_rp microDecayRate(rate_t_rp *r);
rate_t_rp proteinDecayRate(rate_t_rp *r);

rate_t_rp miRAssocConst(rate_t_rp *assoc);
rate_t_rp miRDissocConst(rate_t_rp * dissoc);
rate_t_rp TFAssocConst(rate_t_rp *assoc);
rate_t_rp TFDissocConst(rate_t_rp * dissoc);


unsigned char connsByOutdegree(species_t sp);

UnsignedIntArray getRandomConnectionNumbersForElements(unsigned int leftLen, unsigned int rightLen, species_t left,species_t right);
void initializeParameters(const char * cwd);


#endif