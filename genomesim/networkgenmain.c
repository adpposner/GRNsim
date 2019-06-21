#include "include/modelgenerator.h"
#include "include/modelloader.h"
#include "include/sim.h"
#include "include/globals.h"
#include "include/iofuncs.h"
#include <string.h>
#include <unistd.h>
#include <getopt.h>






int main(int argc, char * argv[]){
	SimulationComponents sim = {0};
    createNetwork(&sim);
	//printDegrees(stdout,sim.g);
	destroyNetwork(&sim);
	return 0;
}
