#include "../include/models/reaction.h"
#include "../include/models/minheap.h"
#include "../include/globals.h"
#include <math.h>
#include <string.h>
// void swapReactionIndices(Reaction * a, Reaction * b) {
// 	unsigned tmp;
// 	tmp = a->index;
// 	a->index = b->index;
// 	b->index = tmp;
// }

void minHeap_free(minHeap * mH){
    if(mH){
        free(mH->indexArr);
        free(mH->rxnArr);
        free(mH->timeToEventArr);
        mH->timeToEventArr = NULL;
        mH->rxnArr = NULL;
        mH->indexArr = NULL;
        mH->len = 0;
    }
}

//DEFINEBASICARRAYTYPEPRIMITIVE(minHeapNode,minHeapNode)
typedef struct minHeapIters{
    ulong_type * const indexStart, *const indexEnd, *indexCurr;
    double * const tteStart, *const tteEnd, *tteCurr;
    Reaction ** const rxnStart, **const rxnEnd, ** rxnCurr;
} minHeapIters;

minHeapIters getMinHeapIters(minHeap * mH){
    return (minHeapIters){.indexStart = mH->indexArr,.indexEnd=mH->indexArr + mH->len,
    .indexCurr = mH->indexArr, .tteStart = mH->timeToEventArr, .tteEnd = mH->timeToEventArr + mH->len,
    .tteCurr = mH->timeToEventArr,.rxnStart = mH->rxnArr, .rxnEnd = mH->rxnArr + mH->len, .rxnCurr = mH->rxnArr};
}

#define MINHEAPITER_FOREACH(it) for(it.indexCurr=it.indexStart,it.rxnCurr=it.rxnStart,it.tteCurr=it.tteStart; \
it.indexCurr != it.indexEnd;it.indexCurr++,it.rxnCurr++,it.tteCurr++) 


extern inline unsigned lesserNodeByIndex(unsigned a, unsigned b, minHeap *mH);
extern inline void swapNodes(minHeap * m, unsigned a, unsigned b);
extern inline unsigned smallestNonParent(minHeap *mH,unsigned p);
extern inline unsigned char pushUp(minHeap * mH,unsigned p);
extern inline unsigned char pushDown(minHeap *mH,unsigned p);






void minHeapify(minHeap *mH, unsigned p) {
	unsigned smallest;
	if (LEFT(p) < mH->len)
	{
		smallest = lesserNodeByIndex(p, LEFT(p), mH);
	}
	else
		smallest = p;
	if (RIGHT(p) < mH->len)
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


static  unsigned char nodeIsFinite(unsigned pos, minHeap * mH)
{return isfinite(mH->timeToEventArr[pos]);}

static  unsigned smallestFiniteChild(minHeap * mH, unsigned p) {

	if (RIGHT(p) < mH->len) {
		if (nodeIsFinite(LEFT(p), mH)) {
			double tL = mH->timeToEventArr[p];
			if (tL < mH->timeToEventArr[RIGHT(p)]) {

				return LEFT(p);
			} else return RIGHT(p);
		} else if (nodeIsFinite(RIGHT(p), mH)) {
			return RIGHT(p);
		}
	} else if (LEFT(p) < mH->len) {
		if (nodeIsFinite(LEFT(p), mH))
			return LEFT(p);
	}
	return p;
}

 unsigned char pushDown_inf(minHeap *mH,unsigned p){
	unsigned curr = p;
	unsigned smallest = curr;
	while(curr!=(smallest = smallestFiniteChild(mH, curr))){
		swapNodes(mH, curr, smallest);
		curr=smallest;
	}
	return curr==p;
}



//  void updateSmall(minHeap * mH, unsigned p) {
// 	assert(mH->data[p].index == p);
// 	unsigned smallest = p;
// 	if (nodeIsFinite(p, mH)) {
// 		if (p != 0) {
// 			if (mH->data[p].timeToEvent < mH->data[PARENT(p)].timeToEvent) {
// 				swapNodes(mH, p, PARENT(p));
// 				updateSmall(mH, PARENT(p));
// 			}
// 		}
// 		smallest = smallestNonParent(mH, smallest);
// 	} else {	//if it's infinite, we just get a smallest finite child
// 		smallest = smallestFiniteChild(mH, smallest);
// 	}
// 	if (smallest != p) {
// 		swapNodes(mH, p, smallest);
// 		updateSmall(mH, smallest);
// 	}
// }



#define PARENTNODE(p)	(p-(p->index/2 +1))
#define LEFTNODE(p)	(p+p->index +1)
#define RIGHTNODE(p) (p+p->index + 2)







 void update(minHeap *mH,unsigned p){
	//updateSmall(mH, p);
	updateSingleNode(mH,p);
}



extern void updateRoot(minHeap *mH);
// 	pushDown(mH, 0);
// }



//Again from Cormen, building the minheap
void buildMinHeap(minHeap *mH) {
	int szHalf;
	unsigned curr;
	szHalf = mH->len-1;//(mH->len / 2) - 1;

	for (curr = szHalf; szHalf >= 0; szHalf--, curr--)
	{
		minHeapify(mH, curr);
	}
}



#include "../include/sim/reactionrates.h"
static void minHeap_index(minHeap *toInit, const ReactionArray *rxns) {
	ReactionArrayIters rIt = getReactionArrayIters(rxns);
	minHeapIters mHIt = getMinHeapIters(toInit);
	unsigned currInd;
    for(mHIt.indexCurr=mHIt.indexStart,mHIt.rxnCurr=mHIt.rxnStart,mHIt.tteCurr=mHIt.tteStart,rIt.curr=rIt.start,currInd=0;
    mHIt.indexCurr != mHIt.indexEnd;mHIt.indexCurr++,mHIt.tteCurr++,mHIt.rxnCurr++,currInd++,rIt.curr++){
        *mHIt.indexCurr = currInd;
        *mHIt.rxnCurr = rIt.curr;
        *mHIt.tteCurr = 0.0;
        
        rIt.curr->timeToEvent = mHIt.tteCurr;
        rIt.curr->mHIndex = *mHIt.indexCurr;
    }
    // }
	// for (rIt.curr = rIt.start, mIt.curr = mIt.start, currInd = 0; mIt.curr != mIt.end; mIt.curr++, rIt.curr++, currInd++) {
	// 	rIt.curr->tte = mIt.curr;
	// 	mIt.curr->index = currInd;
	// 	mIt.curr->rxn = rIt.curr;
	// 	//mIt.curr->timeToEvent = 0.0;
	// }
	//initAllRatesTimes(rxns,toInit);
	assert(rIt.curr == rIt.end);
}

minHeap * minHeap_init(const ReactionArray *rxns) {
	minHeap *toRet;
	toRet = malloc(sizeof(*toRet));
    toRet->len = rxns->length;
    toRet->rxnArr = malloc(toRet->len * sizeof(*toRet->rxnArr));
    toRet->timeToEventArr = malloc(toRet->len * sizeof(*toRet->timeToEventArr));
    toRet->indexArr = malloc(toRet->len * sizeof(*toRet->indexArr));
	minHeap_index(toRet, rxns);
	buildMinHeap(toRet);
	return toRet;
}




//Doesn't pretty print, must fix
void printMinHeap(minHeap mH) {
	minHeapIters it = getMinHeapIters(&mH);
	MINHEAPITER_FOREACH(it){
        printf("("ltf_rp",%u,%f)", *it.indexCurr, (*it.rxnCurr)->id, *it.tteCurr);
    }
		

}

//Helper function to check minHeapProperty
int testMinHeap(minHeap mH) {
	int ctr = 0;
	static int ctrOld = 0;

	unsigned i;
	for (i = 0; i < mH.len; i++)
	{
		assert(i == mH.indexArr[i]);
		if (LEFT(i) < mH.len) {
			//if(mH.data[i]->timeToEvent > mH.data[LEFT(i)]->timeToEvent)
			//{	update(&mH,i);}//update(&mH,LEFT(i));}
			if (mH.timeToEventArr[i] > mH.timeToEventArr[LEFT(i)])
			{	printf("VIOLATIONL\t%d\t%f\t%f\n", i, mH.timeToEventArr[i],
				       mH.timeToEventArr[LEFT(i)]);
				printRxn(stdout, mH.rxnArr[i], 0);
				printRxn(stdout, mH.rxnArr[LEFT(i)], 0);
				ctr++;
			}

		}
		if (RIGHT(i) < mH.len)
		{

			if (mH.timeToEventArr[i] > mH.timeToEventArr[RIGHT(i)])
			{	printf("VIOLATIONR\t%d\t%f\t%f\n", i, mH.timeToEventArr[i],
				       mH.timeToEventArr[RIGHT(i)]);
				printRxn(stdout, mH.rxnArr[i], 0);
				printRxn(stderr, mH.rxnArr[RIGHT(i)], 0); ctr++;
			}
		}
	}
	if (ctr != 0) {
		printf("failure #%d\n", ctrOld); ctrOld++;
	}
	return ctr;

}
