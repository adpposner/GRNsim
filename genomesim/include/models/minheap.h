//minheap.h - Minheap for Gibson-Bruck simulation - modeled using array, not "heap"
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
#ifndef MINHEAP_H_RP__
#define MINHEAP_H_RP__
#include <stdlib.h>
#include "reaction.h"

#define LEFT(i)	(2*i)+1
#define RIGHT(i) 2*(i+1)
#define PARENT(i) (i-1)/2


typedef struct minHeapNode {
	ulong_type index;
	double timeToEvent;
	Reaction * rxn;
} minHeapNode;

DECLAREBASICARRAYTYPEPRIMITIVE(minHeapNode,minHeapNode)

typedef minHeapNodeArray minHeap;

minHeap *minHeap_init(const ReactionArray *rxns);
void minHeap_free(minHeap *mH);
void minHeapify(minHeap * mH,unsigned p);
void buildMinHeap(minHeap * mH);

void updateSmall(minHeap * mH, unsigned p);
void update(minHeap * mH,unsigned p);

//
inline unsigned lesserNode(minHeapNode *a, minHeapNode * b) {

	return a->timeToEvent > b->timeToEvent ? b->index : a->index;
}

//
inline unsigned lesserNodeByIndex(unsigned a, unsigned b, minHeap *mH) {
	return lesserNode(mH->data + a, mH->data + b);
}


//
inline void swapNodes(minHeap * m, unsigned a, unsigned b) {
 	minHeapNode temp = m->data[a];
 	m->data[a]=m->data[b];
 	m->data[b]=temp;
 	m->data[a].index = a;
 	m->data[b].index = b;
 	m->data[a].rxn->tte = &m->data[a];
 	m->data[b].rxn->tte = &m->data[b];
 }

//
inline unsigned smallestNonParent(minHeap *mH,unsigned p){
	unsigned smallest = p;
	if (RIGHT(p) < mH->length) {
		smallest = lesserNodeByIndex(p, LEFT(p), mH);
		smallest = lesserNodeByIndex(smallest, RIGHT(p), mH);
	} else if (LEFT(p) < mH->length) {
		smallest = lesserNodeByIndex(smallest, LEFT(p), mH);
	}
	return smallest;
}

// Double declaration for C99 std inlining rules
inline unsigned char pushUp(minHeap * mH,unsigned p){
	unsigned curr = p;
	while (curr && mH->data[curr].timeToEvent < mH->data[PARENT(curr)].timeToEvent){
		swapNodes(mH, curr, PARENT(curr));
		curr=PARENT(curr);
	}
	return curr == p;
}

// Double declaration for C99 std inlining rules
inline unsigned char pushDown(minHeap *mH,unsigned p){
	unsigned curr = p;
	unsigned smallest = curr;
	while (curr!=(smallest = smallestNonParent(mH, curr))){
		swapNodes(mH, curr, smallest);
		curr=smallest;
	}	
	return curr == p;
}



// unsigned char pushUp(minHeap * mH,unsigned p);
// unsigned char pushDown(minHeap *mH,unsigned p);

//
inline void updateRoot(minHeap * mH){	
	pushDown(mH, 0);
}

void printMinHeap(minHeap mH);
int testMinHeap(minHeap mH);




#endif