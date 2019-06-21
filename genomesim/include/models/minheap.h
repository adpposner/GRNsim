#ifndef __MINHEAP_H__
#define __MINHEAP_H__
#include <stdlib.h>
#include "reaction.h"

#define LEFT(i)	(2*i)+1
#define RIGHT(i) 2*(i+1)
#define PARENT(i) (i-1)/2


// typedef struct minHeapNode {
// 	ulong_type index;
// 	double timeToEvent;
// 	Reaction * rxn;
// } minHeapNode;

typedef struct minHeapSoA {
    ulong_type * indexArr;
    double * timeToEventArr;
    Reaction ** rxnArr;
    ulong_type len;
} minHeapSoA;

//DECLAREBASICARRAYTYPEPRIMITIVE(minHeapNode,minHeapNode)

typedef minHeapSoA minHeap;

minHeap *minHeap_init(const ReactionArray *rxns);
void minHeap_free(minHeap *mH);
void minHeapify(minHeap * mH,unsigned p);
void buildMinHeap(minHeap * mH);

//void updateSmall(minHeap * mH, unsigned p);
void update(minHeap * mH,unsigned p);

// inline unsigned lesserNode(minHeapNode *a, minHeapNode * b) {
// 	// double at = a->timeToEvent, bt = b->timeToEvent;
// 	return a->timeToEvent > b->timeToEvent ? b->index : a->index;
// }

inline unsigned lesserNodeByIndex(unsigned a, unsigned b, minHeap *mH) {
    return (mH->timeToEventArr[a] <= mH->timeToEventArr[b]) ? mH->indexArr[a] : mH->indexArr[b];
}



inline void swapNodes(minHeap * m, unsigned a, unsigned b) {
    ulong_type indexTemp;
    Reaction * rxnTemp;
    double tteTemp;
    assert(a == m->indexArr[a]);
    assert(b == m->indexArr[b]);
    assert(a == m->rxnArr[a]->mHIndex);
    assert(b == m->rxnArr[b]->mHIndex);
    tteTemp = m->timeToEventArr[a];
    rxnTemp = m->rxnArr[a];
    indexTemp = m->indexArr[a];
    //m->indexArr[a] = m->indexArr[b];
    m->timeToEventArr[a] = m->timeToEventArr[b];
    m->rxnArr[a] = m->rxnArr[b];
    m->rxnArr[a]->timeToEvent = &m->timeToEventArr[a];
    m->rxnArr[a]->mHIndex = m->indexArr[a];
    //m->indexArr[b] = indexTemp;
    m->timeToEventArr[b] = tteTemp;
    m->rxnArr[b] = rxnTemp;
    m->rxnArr[b]->timeToEvent = &m->timeToEventArr[b];
    m->rxnArr[b]->mHIndex = m->indexArr[b];
 }


inline unsigned smallestNonParent(minHeap *mH,unsigned p){
	unsigned smallest = p;
	if (RIGHT(p) < mH->len) {
		smallest = lesserNodeByIndex(p, LEFT(p), mH);
		smallest = lesserNodeByIndex(smallest, RIGHT(p), mH);
	} else if (LEFT(p) < mH->len) {
		smallest = lesserNodeByIndex(smallest, LEFT(p), mH);
	}
	return smallest;
}

inline unsigned char pushUp(minHeap * mH,unsigned p){
	unsigned curr = p;
    while(curr &&( mH->timeToEventArr[curr] < mH->timeToEventArr[PARENT(curr)])){
        swapNodes(mH,curr,PARENT(curr));
        curr = PARENT(curr);
    }
	return curr == p;
}

inline unsigned char pushDown(minHeap *mH,unsigned p){
	unsigned curr = p;
	unsigned smallest = curr;
	while (curr!=(smallest = smallestNonParent(mH, curr))){
		swapNodes(mH, curr, smallest);
		curr=smallest;
	}	
	return curr == p;
}



unsigned char pushUp(minHeap * mH,unsigned p);
unsigned char pushDown(minHeap *mH,unsigned p);

//void updateRoot(minHeap *mH);
inline void updateRoot(minHeap * mH){	
	pushDown(mH, 0);
}

void printMinHeap(minHeap mH);
int testMinHeap(minHeap mH);

void generate_from_mh(minHeap mH,Reaction * toMove);




#endif