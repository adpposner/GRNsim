// minheap.c - contains all minheap update funcs
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
#include "../include/models/reaction.h"
#include "../include/models/minheap.h"
#include "../include/globals.h"
#include <math.h>
#include <string.h>

void minHeapNode_free(minHeapNode * mH) {}

DEFINEBASICARRAYTYPEPRIMITIVE(minHeapNode,minHeapNode)

//Extern declarations for C99 inlining std
extern inline unsigned lesserNode(minHeapNode *a, minHeapNode * b);
extern inline unsigned lesserNodeByIndex(unsigned a, unsigned b, minHeap *mH);
extern inline void swapNodes(minHeap * m, unsigned a, unsigned b);
extern inline unsigned smallestNonParent(minHeap *mH,unsigned p);
extern inline unsigned char pushUp(minHeap * mH,unsigned p);
extern inline unsigned char pushDown(minHeap *mH,unsigned p);



void minHeapify(minHeap *mH, unsigned p) {
	unsigned smallest;
	if (LEFT(p) < mH->length)
	{
		smallest = lesserNodeByIndex(p, LEFT(p), mH);
	}
	else
		smallest = p;
	if (RIGHT(p) < mH->length)
		smallest = lesserNodeByIndex(smallest, RIGHT(p), mH);
	if (smallest != p)
	{
		swapNodes(mH, p, smallest);
		minHeapify(mH, smallest);
	}
}




 void updateSingleNode(minHeap *mH,unsigned p){
	unsigned didntPushUp = pushUp(mH, p);
	if (didntPushUp){
		pushDown(mH, p);
	}
}



 void update(minHeap *mH,unsigned p){

	updateSingleNode(mH,p);
}



extern void updateRoot(minHeap *mH);




void buildMinHeap(minHeap *mH) {
	int szHalf;
	unsigned curr;
	szHalf = mH->length-1;

	for (curr = szHalf; szHalf >= 0; szHalf--, curr--)
	{
		minHeapify(mH, curr);
	}
}



#include "../include/sim/reactionrates.h"
static void minHeap_index(minHeap *toInit, const ReactionArray *rxns) {
	ReactionArrayIters rIt = getReactionArrayIters(rxns);
	minHeapNodeArrayIters mIt = getminHeapNodeArrayIters(toInit);
	unsigned currInd;
	for (rIt.curr = rIt.start, mIt.curr = mIt.start, currInd = 0; mIt.curr != mIt.end; mIt.curr++, rIt.curr++, currInd++) {
		rIt.curr->tte = mIt.curr;
		mIt.curr->index = currInd;
		mIt.curr->rxn = rIt.curr;

	}

	assert(rIt.curr == rIt.end);
}

minHeap * minHeap_init(const ReactionArray *rxns) {
	minHeap *toRet;
	toRet = malloc(sizeof(*toRet));
	*toRet = minHeapNodeArray_alloc(rxns->length);
	minHeap_index(toRet, rxns);
	buildMinHeap(toRet);
	return toRet;
}


void minHeap_free(minHeap *mH) {
	if(mH)
		minHeapNodeArray_free(mH);
	mH=NULL;
}

//Bare printing of heap - not a "pretty print" functions
void printMinHeap(minHeap mH) {
	minHeapNodeArrayIters it = getminHeapNodeArrayIters(&mH);
	for (it.curr = it.start; it.curr != it.end; it.curr++)
		printf("("ltf_rp",%u,%f)", it.curr->index, it.curr->rxn->id, it.curr->timeToEvent);

}

//Helper function to check minHeapProperty
int testMinHeap(minHeap mH) {
	int ctr = 0;
	static int ctrOld = 0;

	unsigned i;
	for (i = 0; i < mH.length; i++)
	{
		assert(i == mH.data[i].index);
		if (LEFT(i) < mH.length) {

			if (mH.data[i].timeToEvent > mH.data[LEFT(i)].timeToEvent)
			{	printf("VIOLATIONL\t%d\t%f\t%f\n", i, mH.data[i].timeToEvent,
				       mH.data[LEFT(i)].timeToEvent);
				printRxn(stdout, mH.data[i].rxn, 0);
				printRxn(stdout, mH.data[LEFT(i)].rxn, 0);
				ctr++;
			}

		}
		if (RIGHT(i) < mH.length)
		{

			if (mH.data[i].timeToEvent > mH.data[RIGHT(i)].timeToEvent)
			{	printf("VIOLATIONR\t%d\t%f\t%f\n", i, mH.data[i].timeToEvent,
				       mH.data[RIGHT(i)].timeToEvent);
				printRxn(stdout, mH.data[i].rxn, 0);
				printRxn(stderr, mH.data[RIGHT(i)].rxn, 0); ctr++;
			}
		}
	}
	if (ctr != 0) {
		printf("failure #%d\n", ctrOld); ctrOld++;
	}
	return ctr;

}
