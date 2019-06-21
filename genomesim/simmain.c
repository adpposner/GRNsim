#include "include/modelloader.h"
#include "include/sim.h"
#include "include/globals.h"
#include <string.h>
#include <unistd.h>
#include <getopt.h>





SimOpts parseOpts(int argc, char * argv[]){
	SimOpts toRet = {.jsonFile_ = '\0', .dataOutputDir_ = '\0',.optSuffix_ = '\0',
     .writeGenes = 1, .nactiveMirs=-1,.nSimulations =-1,};

	int opt;
	while ((opt = getopt(argc, argv, "f:gn:s:o:"))!=-1){
		switch(opt){
			case 'f': 	strncpy(toRet.jsonFile_, optarg, MAX_JSON_FILENAME);
						break;
			case 'n':	toRet.nactiveMirs = atoi(optarg);break;
			case 's':	toRet.nSimulations = atoi(optarg);break;
            case 'o':   strncpy(toRet.dataOutputDir_,optarg,MAX_JSON_FILENAME);break;
			default: fprintf(stderr,"Usage %s [-f jsonFileInput - required] [-o destinationDirectory] [-s nSimulations] [-n nActiveMirs] \n",argv[0]);
			exit(EXIT_FAILURE);
		}
	}
    if (!strlen(toRet.jsonFile_)){
        fprintf(stderr,"No filename specified, use -f [jsonFileInput] \n");
        exit(EXIT_FAILURE);
    }
    if (!strlen(toRet.dataOutputDir_)){
        fprintf(stderr,"No output directory specified, defaulting to , genes list may overwrite contents\n");
        snprintf(toRet.dataOutputDir_,MAX_JSON_FILENAME,"output");
    }
    if (toRet.nactiveMirs <= 0){
        fprintf(stderr,"No number of active miRNA specified, defaulting to 0\n");
        toRet.nactiveMirs = 0;
    }
    if (toRet.nSimulations <= 0){
        fprintf(stderr,"No number of simulations specified, defaulting to 1\n");
        toRet.nSimulations = 1;
    }
	return toRet;
}





int main(int argc, char * argv[]){
	SimulationComponents sim = {0};
	SimOpts sO = parseOpts(argc, argv);
 
	loadModel(&sim, &sO);
	doSimulation(&sim,sO.nSimulations);
	destroyModel(&sim);
	return 0;
}
