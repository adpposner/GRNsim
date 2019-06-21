# GRNsim
A high-performance stochastic simulator and network generator for discrete gene regulatory network model with TF and miRNA-mediated gene regulation

Note: This document and source are currently intended for reviewers of the associated miRNA-control manuscript/submission by R. Posner and R. Laubenbacher. This is not a "general-purpose" readme at the moment, but will be expanded accordingly following review and/or acceptance.

Three programs are contained in this repository. They are:

* **Network generator** - used for creation of GRNs according to a set of specified parameters. See `netgen.md`
* **Network simulator** - core program. Simulates regulated transcription and translation for any (not just specified by network generator) GRN in JSON format. See `netinput.md`
* **Simulation analysis** - tool written in C and Fortran for calculating times to differentiation and other dynamic/statistical properties. See `simanalysis/README.md`



## Acknowledgments

* The JSON file reading subroutines (genomesim/include/cjson.h and genomesim/iofuncs/cjson.c) were written by Dave Gamble and are from his (excellent) cJSON project at https://github.com/DaveGamble/cJSON
* While this has been tested (and works) using the gcc's C compiler and clang, it has been optimized for use with the Intel C and Intel Fortran Compilers (2016-2018 versions, https://software.intel.com/en-us/intel-compilers). 
* This project currently uses the Intel Math Kernel Library (https://software.intel.com/mkl) for pseudorandom number generation, FFTs, and BLAS operations. 
