#ifndef GENESLIST_P_H__
#define GENESLIST_P_H__

#include "basicgeneticelement.h"

typedef struct initialQuantities{
	UnsignedIntArray prodICs;
	UnsignedIntArray modICs;
	UnsignedIntArray boundICs;
} initialQuantities;

typedef struct GenesList {
	ProducerArray producers;
	ModulatorArray modulators;
	BoundElementArray bounds;
	ulong_type nDNA, nMess, nMicro,nProt,nMessMir,nTFDNA;
	initialQuantities * ICs;
} GenesList;


typedef struct GenesListIterator {
	ProducerArrayIters pIt;
	ModulatorArrayIters mIt;
	BoundElementArrayIters bIt;
} GenesListIterator;

GenesListIterator getGenesListIters(GenesList * g);
void extractInitialQuantities(GenesList * g);

#endif