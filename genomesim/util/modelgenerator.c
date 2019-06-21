#include "../include/modelgenerator.h"
#include "../include/globals.h"
#include "../include/randomnumbers.h"
#include "../include/parameters.h"
#include "../include/sim.h"
#include "../include/iofuncs.h"
#include "../include/models/geneslist_p.h"
#include "../include/models/basicgeneticelement_p.h"
#include "../include/models/reaction.h"
#include <string.h>









//want to package below so the indegree/outdegree things aren't split
static unsigned extendBoundElementArray(){return 0;}





static ulong_type createRandomConnections_(GenesList *g, species_t left, species_t right) {
	if(left & DNA) assert(right & PROTEIN);
	if(left & MESSENGER) assert(right & MICRO);
	pmbPtrArray leftPtrArr;
	pmbPtrArray rightPtrArr;
	leftPtrArr = getPtrArrayForType(g, left);
	rightPtrArr = getPtrArrayForType(g, right);


	UnsignedIntArray nConns;
	int by_indegree=0;
	//Get the indegrees for each element - !!Using indegrees

	nConns = getRandomConnectionNumbersForElements(leftPtrArr.length,rightPtrArr.length, left,right);

	unsigned totalConns;
	//Sum up to find out how far we have to extend our arrays
	totalConns = UnsignedIntArray_total(nConns);
	//are we just adding more boundelements?
	ulong_type oldLength = g->bounds.length;
	BoundElementArray_extend(&g->bounds, totalConns);
	//If we've already done our first set of connections (i.e. if we are now doing miRs vs. )
	BoundElementArrayIters bIt = getBoundElementArrayIters(&g->bounds);
	bIt.curr = bIt.start + oldLength;
	
	pmbPtrArray *source = NULL,*target = NULL;
	
	if(connsByOutdegree(left) == INDEGREE){
		source = &leftPtrArr;
		target = &rightPtrArr;

	pmbPtrArrayIters srcIt = getpmbPtrArrayIters(source);
	pmbPtrArrayIters targIt = getpmbPtrArrayIters(target);
	UnsignedIntArrayIters connIt = getUnsignedIntArrayIters(&nConns);
	
	unsigned k;
	connIt.curr=connIt.start;
	ARRAY_TYPE_FOREACH(srcIt){

		shufflePtrArray(target->data, target->length, sizeof(*target->data));
		targIt.curr = targIt.start;
		for(k=0;k<*connIt.curr;k++){
			
			bIt.curr = initConnection(srcIt.curr, targIt.curr++, bIt.curr);
		}
		connIt.curr++;
	}
	}else if(connsByOutdegree(right) == OUTDEGREE){
		source = &rightPtrArr;
		target = &leftPtrArr;

	pmbPtrArrayIters srcIt = getpmbPtrArrayIters(source);
	pmbPtrArrayIters targIt = getpmbPtrArrayIters(target);
	UnsignedIntArrayIters connIt = getUnsignedIntArrayIters(&nConns);
	
	unsigned k;
	connIt.curr=connIt.start;
	ARRAY_TYPE_FOREACH(srcIt){

		shufflePtrArray(target->data, target->length, sizeof(*target->data));
		targIt.curr = targIt.start;
		for(k=0;k<*connIt.curr;k++){
			
			bIt.curr = initConnection(targIt.curr++, srcIt.curr, bIt.curr);
		}
		connIt.curr++;
	}
	}
	assert(bIt.curr==bIt.end);
	pmbPtrArray_free(&leftPtrArr);
	pmbPtrArray_free(&rightPtrArr);
	UnsignedIntArray_free(&nConns);
	return totalConns;
}


static void randomlyConnectGeneslist(GenesList *g) {
	unsigned nMessMir = 0, nTFDNA = 0;
	nMessMir = createRandomConnections_(g,MESSENGER,MICRO);
	nTFDNA = createRandomConnections_(g,CODING, PROTEIN);
	nTFDNA += createRandomConnections_(g,NONCODING,PROTEIN);
	genesList_bound_quantities_init(g, nMessMir, nTFDNA);
	assert((g->nMessMir + g->nTFDNA) == g->bounds.length);
	assignBoundElements(g);
}

GenesList * setupRandomGenesList(){

	GenesList *toRet = genesList_base_generate(globalDims.nMess,globalDims.nMicro);
	randomlyConnectGeneslist(toRet);
	return toRet;
}


void createNetwork(SimulationComponents *sim){
    	initializeParameters(".");
		randomnumbers_init(networkGenSeed);
		sim->g = setupRandomGenesList();
		initDefaultGenesListQuantities(sim->g);
		writeGenesList(sim->g,INCLUDE_DISABLED,NULL,0,NULL);
}

void destroyNetwork(SimulationComponents *toFree) {
    destroyParameters();
	minHeap_free(toFree->mH);
	if(toFree->rxns)
		ReactionArray_free(toFree->rxns);
	toFree->rxns = NULL;
	if(toFree->g)
		genesList_free(toFree->g);
	toFree->g = NULL;
	randomnumbers_free();
}