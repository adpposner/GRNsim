#include "../include/iofuncs/iofuncs.h"
#include "../include/iofuncs/graph_interface.h"
#include "../include/models/basicgeneticelement.h"
#include "../include/models/basicgeneticelement_p.h"
#include <cgraph.h>
#include <gvc.h>



//Agraph_t *G;
//
//void generate_graph() {
//	G = agopen("network", Agstrictdirected, 0);
//	agsafeset(G, "nodesep", "0.5", "1.0");
//	agsafeset(G, "rankdir", "LR", "LR");
//
//}
//
//// Agnode_t * generateNodeForElem(BasicGeneticElement * bge) {
//// 	char idbuf[100];
//// 	Agnode_t *toRet;
//// 	sprintf(idbuf,"%u",bge->id);
//
//// 	Agraph_t * gene;
//// 	char clusterName[100];
//// 	sprintf(clusterName,"cluster%s",bge->name);
//// 	gene = agsubg(G,clusterName,1);
//// 	toRet = agnode(gene, idbuf,1);
//// 	agsafeset(toRet,"label","","");
//// 	agsafeset(gene,"rankdir","TD","TD");
//// 	agsafeset(gene,"rank","same","same");
//
//// 	agsafeset(gene,"label",bge->name,"");
//// 	//agsafeset(gene,"rankdir","LR","LR");
//// 	if (bge->species & (MICRO | MESSENGER))
//// 		agsafeset(toRet,"shape","triangle","ellipse");
//// 	else if (bge->species & (PROTEIN))
//// 		agsafeset(toRet,"shape","square","square");
//// 		//agsafeset(toRet,"shape","square","ellipse");
//// 	return toRet;
//// }
//
//Agnode_t * generateShortHandFor(pmbPtr p) {
//	char idbuf[100];
//	Agnode_t *toRet = NULL;
//
//	if (p.species & (DNA | MESSENGER)) {
//		sprintf(idbuf, "%zu", p.p.producer->info.id);
//		toRet = agnode(G, idbuf, 1);
//
//		agsafeset(toRet, "label", p.p.producer->info.name, "");
//	} else if (p.species & (PROTEIN | MICRO)) {
//		sprintf(idbuf, "%zu", p.p.modulator->info.id);
//		toRet = agnode(G, idbuf, 1);
//		agsafeset(toRet, "label", p.p.modulator->info.name, "");
//		if (p.species & (MICRO)) 
//			agsafeset(toRet, "shape", "ellipse", "ellipse");
//		else if (p.species & (PROTEIN)) {
//			agsafeset(toRet, "shape", "square", "square");
//		}		//agsafeset(toRet,"shape","square","ellipse");
//}
//return toRet;
//}
//
//
//void generateShorthandEdges(BasicGeneticElementArray bd) {
//	BasicGeneticElementArrayIters it = getBasicGeneticElementArrayIters(&bd);
//	char edgename[100];
//	char headbuf[100];
//	char tailbuf[100];
//	Agnode_t * head = NULL, *tail = NULL;
//	Agedge_t * edge;
//	for (it.curr = it.start; it.curr != it.end; it.curr++) {
//		sprintf(tailbuf, "%s", it.curr->base.bound.left.elem->name);
//
//		sprintf(headbuf, "%s", it.curr->base.bound.right.elem->name);
//		sprintf(edgename, "%u", it.curr->id);
//		head = agnode(G, headbuf, 0);
//		tail = agnode(G, tailbuf, 0);
//		edge = agedge(G, head, tail, edgename, 1);
//		if (it.curr->base.bound.effectStrength < 1.0)
//			agsafeset(edge, "arrowhead", "tee", "tee");
//		else agsafeset(edge, "arrowhead", "normal", "normal");
//	}
//}
//
//void generateProducerNodes(ProducerArray *p) {
//
//	char idB
//	ProducerArrayIters pIt = getProducerArrayIters(p);
//	ARRAY_TYPE_FOREACH(pIt) {
//
//	}
//}
//
//void generateBoundNodes(BasicGeneticElementArray bd) {
//	BasicGeneticElementArrayIters it = getBasicGeneticElementArrayIters(&bd);
//	char edgename[100];
//	char headbuf[100];
//	char tailbuf[100];
//	Agnode_t * head, *tail;
//	Agedge_t * edge;
//	for (it.curr = it.start; it.curr != it.end; it.curr++) {
//		sprintf(tailbuf, "%u", it.curr->base.bound.left.elem->id);
//		sprintf(headbuf, "%u", it.curr->base.bound.right.elem->id);
//		sprintf(edgename, "%u", it.curr->id);
//		head = agnode(G, headbuf, 0);
//		tail = agnode(G, tailbuf, 0);
//		edge = agedge(G, head, tail, edgename, 1);
//		if (it.curr->base.bound.effectStrength < 1.0)
//			agsafeset(edge, "arrowhead", "tee", "tee");
//		else agsafeset(edge, "arrowhead", "normal", "normal");
//	}
//}
//
//void generateEdgesForElem(BasicGeneticElement * bge) {
//	char headbuf[100], tailbuf[100];
//	if (isProducer(bge)) {
//		sprintf(headbuf, "%u", bge->id);
//		sprintf(tailbuf, "%u", bge->base.base.produces.elem->id);
//		Agnode_t * head, *tail;
//		Agedge_t *edge;
//		head = agnode(G, headbuf, 0);
//		tail = agnode(G, tailbuf, 0);
//
//		edge = agedge(G, head, tail, NULL, 1);
//		agsafeset(edge, "arrowhead", "normal", "normal");
//	}
//}
//
//
//void graph_toFile(const char * filename) {
//
//	GVC_t * gvc;
//	gvc = gvContext();
//	agsafeset(G, "overlap", "false", "false");
//
//	gvLayout(gvc, G, "dot");
//
//	gvRenderFilename(gvc, G, "png", filename);
//	gvFreeLayout(gvc, G);
//}
//
//void shorthandGraphForBGEArray(connectedGenesList *cg, const char * filename) {
//	BasicGeneticElementArrayIters it = getBasicGeneticElementArrayIters(&cg->genes);
//	generate_graph();
//	for (it.curr = it.start; it.curr != it.end; it.curr++)
//		generateShortHandForElement(it.curr);
//	for (it.curr = it.start; it.curr != it.end; it.curr++)
//	{	generateShorthandEdges(it.curr->base.base.boundelts.elements);}
//
//	graph_toFile(filename);
//	agclose(G);
//}
//
//void graphForBGEArray(connectedGenesList *cg, const char * filename) {
//	BasicGeneticElementArrayIters it = getBasicGeneticElementArrayIters(&cg->genes);
//	generate_graph();
//	for (it.curr = it.start; it.curr != it.end; it.curr++)
//		generateNodeForElem(it.curr);
//	for (it.curr = it.start; it.curr != it.end; it.curr++)
//	{	generateEdgesForElem(it.curr); generateBoundNodes(it.curr->base.base.boundelts.elements);}
//
//	graph_toFile(filename);
//	agclose(G);
//}
