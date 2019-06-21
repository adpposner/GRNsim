#ifndef REACTANT_REF_H__
#define REACTANT_REF_H__

#include "../globals.h"
#include "basicgeneticelement.h"




typedef struct GenesList GenesList;



void resetInitialQuantities(GenesList *g);
pmbPtrArray getPtrArrayForType(GenesList *g, species_t sptype);


GenesList * genesList_alloc(const ulong_type nMess, const ulong_type nMicro, 
	const ulong_type nMessMir, const ulong_type nTFDNA,unsigned char allocInternals);

GenesList * genesList_base_generate(ulong_type nMess, ulong_type nMicro);
void genesList_bound_quantities_init(GenesList *g, const ulong_type nMessMir, const ulong_type nTFDNA);
void assignBoundElements(GenesList *g);
void initDefaultGenesListQuantities(GenesList *g);

void genesList_free(GenesList * tofree);

ulong_type nMessengers(GenesList * g);

ulong_type nMicro_total(GenesList * g);

ulong_type nEnabledForType(GenesList * g,species_t species);

void toggleFamilyForDisabled(GenesList *g);
unsigned isEnabled(const pmbPtr pmb);
unsigned char isProducerEnabled(Producer * producer);
unsigned char isModulatorEnabled(Modulator * modulator);
unsigned char isBoundEnabled(BoundElement * bd);
#endif