// mainnetgen.c - main function for network generator
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
#include "include/modelgenerator.h"
#include "include/modelloader.h"
#include "include/sim.h"
#include "include/globals.h"
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#include <getopt.h>





SimOpts parseOpts(int argc, char * argv[]){
	SimOpts toRet = {0};
    *toRet.jsonFile_ = '\0';

	int opt;
	while ((opt = getopt(argc, argv, "o:"))!=-1){
		switch(opt){
			case 'o': 	strncpy(toRet.jsonFile_, optarg, MAX_JSON_FILENAME);
						break;
			default: fprintf(stderr,"Usage %s -o [destinationJSONpath] \n",argv[0]);
			exit(EXIT_FAILURE);
		}
	}
    if (!strlen(toRet.jsonFile_)){
        fprintf(stderr,"No filename specified, use -o [OUTFILENAME] \n");
        exit(EXIT_FAILURE);
    }
	return toRet;
}


int main(int argc, char * argv[]){
	SimulationComponents sim = {0};
    SimOpts opts = {0};
    opts = parseOpts(argc,argv);
    createNetwork(&sim,opts.jsonFile_);

	destroyNetwork(&sim);
	return 0;
}
