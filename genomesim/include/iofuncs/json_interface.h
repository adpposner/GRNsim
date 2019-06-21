#ifndef JSON_INTERFACE_H__
#define JSON_INTERFACE_H__
#include "cjson.h"
#include "../globals.h"
// cJSON * producer_toJSON(Producer * elem);

// Producer producer_fromJSON(cJSON * obj);
// char * producer_print(Producer * p);



struct GenesList;

void writeGenesList_JSON(struct GenesList * g,const char * filename, SkipDisabled skipEnabled);
struct GenesList * readGenesList_JSON(const char * filename);
//struct genesList readGenesList_JSON(const char * filename);
//void writeGenesList_JSON(struct genesList * g,const char * filename);

#endif