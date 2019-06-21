# GRNsim
A high-performance stochastic simulator and network generator for discrete gene regulatory network model with TF and miRNA-mediated gene regulation

Note: This document and source are currently intended for reviewers of the associated miRNA-control manuscript/submission by R. Posner and R. Laubenbacher. This is not a "general-purpose" readme at the moment, but will be expanded accordingly following review and/or acceptance.

Three programs are contained in this repository. They are:

* **Network generator** - used for creation of GRNs according to a set of specified parameters. See `netgen.md`
* **Network simulator** - core program. Simulates regulated transcription and translation for any (not just specified by network generator) GRN in JSON format. See `netinput.md`
* **Simulation analysis** - tool written in C and Fortran for calculating times to differentiation and other dynamic/statistical properties. See `simanalysis/README.md`

This branch contains a version of the network generator/simulator which has been adapted so that its only dependencies are on free, open-source tools. This version should be compatible with most 64-bit flavors of Linux and Mac OS. It additionally contains a Dockerfile to build a Docker container to build and run the system. One can either build from Docker or build manually. The Fortran/C subroutine source is also available here.

## Building
### Build from Docker
The included Dockerfile can be used to build a container suitable for building the project, and will build the project therein. It should be as simple as entering

>`docker build -t [IMAGE_NAME] .`

The container can then be run using

>`docker run -it [IMAGE_NAME]`

which will open a shell from which the binaries can be run. If one wants to persist the data, instead one can use:

>`docker run -it -v "$(pwd)"/:genomesim [IMAGE_NAME]`

### Build from source (non-Docker)
To build the source, one needs either the GNU C Compiler (GCC) or clang and GNU make. The network generator also depends on *libxml2*, which is installed on Mac by default and is available on all common Linux package managers (apt, yum, pacman, etc).

To build it, ensure that the include path has been set to find libxml2 in `genomesim/util/Makefile` and that `$(CC)` is set accordingly in `genomesim/Makefile`. Beyond that, from the `genomesim` directory, one can simply type:

-`make` to build the system with no miRNA recycling 
-`make recyc` to build the system with miRNA recycling 
-`make clean` to delete build files/binaries 
-`make debug` which includes extra runtime assertipons and tests to ensure validity.

## Usage
### Running the network generator
The network generator has been stripped down to produce single GRNs in JSON format which can be read by the simulator. It pulls its parameters from `genomesim/config.xml`, as specified in the manuscript. Beyond that,

>`./mainnetgen -o [OUTPUTFILENAME]`

will generate a network and place it in the path given in `OUTPUTFILENAME`. See `netgen.md` for an overview of the components of the configuration file. Note that this version is limited to a maximum of 63 binding sites per gene/mRNA - the complete branch has a maximum of 255, as used in the manuscript.

### Running the simulator
The simulator takes the generated JSON file as input. It can be run by entering:

>`./mainsim -f [INPUTJSONFILE] -o [OUTPUTDIRECTORY] -n [N_ACTIVE_MIRNA] -s [NUMBER_OF_RUNS]`

These parameters correspond to the input file path, the directory to output simulation results, the number of miRNAs in the system to activate, and the number of simulation trials. If not specified, the output directory is `output` by default. The simulator will select `N_ACTIVE_MIRNA`s to activate at random, and run that network `NUMBER_OF_RUNS` times. This set of runs will have the exact same set of miRNAs activated - to test different sets, the full simulation command should be invoked again.

### Simulation output
The system will pseudorandomly generate a seed for the simulation, specified below by `SIMNO`. This can be used with in conjunction with the JSON file listed below to reproduce the exact same simulation run. In the specified output directory, one will find four files. These are:

>`INPUTJSONFILE` - from above, the file name is repeated and placed in the output directory. This file corresponds to the original transcriptional network, as well as the n pseudorandomly-selected miRNAs. This will be overwritten if a different output directory is not specified.

>`sim.[SIMNO].prot.txt` - This list contains TF copy numbers at approximately uniformly-spaced time intervals. The header is tab-delimited, with clockTime in the first column, and each protein name subsequently listed.

>`sim.[SIMNO].rna.txt` - This list is similar to the previous, but with mRNA copy numbers.

>`sim.[SIMNO].prot.degrees.txt` - The calculated mRNA indegrees (i.e. how many miRNAs target each protein's mRNA transcript) for the run.

Analysis
The analysis routines are originally written in Fortran with a C interface to achieve the throughput needed for the generation of our results. However, included in the `genomesim` directory one can find a small python script which will do a quick time-to-differentiation evaluation. It is run using:

>`./simpleanalysis.py [PROTEINFILEPATH]`
where `PROTEINFILEPATH` is the emitted `\*.prot.txt` file from the output. This does depend on the numpy and pandas libraries (pandas only for loading the file), which can be easily installed on any python distribution.




## Acknowledgments

* The JSON file reading subroutines (genomesim/include/cjson.h and genomesim/iofuncs/cjson.c) were written by Dave Gamble and are from his (excellent) cJSON project at https://github.com/DaveGamble/cJSON
* While this has been tested (and works) using the gcc's C compiler and clang, it has been optimized for use with the Intel C and Intel Fortran Compilers (2016-2018 versions, https://software.intel.com/en-us/intel-compilers). 
* This project currently uses the Intel Math Kernel Library (https://software.intel.com/mkl) for pseudorandom number generation, FFTs, and BLAS operations. 
