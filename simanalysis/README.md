# Simulation analysis subroutines for GRN simulation

Note: This document and source are currently intended for reviewers of the associated miRNA-control manuscript/submission by R. Posner and R. Laubenbacher. This is not a "general-purpose" readme at the moment, but will be expanded accordingly following review and/or acceptance.

The analysis subroutine and `simanalysis` program are intended for high-throughput processing of GRN simulator data.

## Invocation and input
Call as `main [outputpath] [isProt]` where `outputpath` and `isProt` are required. `outputpath` should be the path of simulation output - i.e. the parent directory of the `net.md.[MD5]` directories created using the network generator, typically just called `[mainsimdirectory]/output`. `isProt` tells the tool whether to calculate time-to-differentiation using discrete Hilbert transform (DHT) of mRNA (when set to 0) or TF (when set to 1) simulation trajectories.

## Output

Output for all processed data is put into a single tab-delimited file called `processedauto.[prot|rna].txt` in the current working directory. The determination of time-to-differentiation depends on whether `isProt` is set. If it is not set, the time-to-differentiation will be calculated using RNA expression levels.

The output file contains the following fields, which are grouped for brevity. values called `xxxx[nElem]` correspond to tab-delimited arrays containing `nElem` items.

* netmd, ngs, nmicro, states, simno, nElem: Grouping variables extracted from file path, corresponding to MD5 parameter hash, network generation seed, total number of active miRNAs, active/inactive states of first 40 miRNAs, simulation seed, and \# of elements (i.e. number of TF species)
* globalttd: Global time to differentiation in hours. Uses DHT of summed, normed copy \# vectors (as in manuscript).
* messPre, microPre, protPre: average mRNA, miRNA, and TF synthesis *rate* (molecules per hour) prior to differentiation.
* messPost, microPost, protPost: average mRNA, miRNA, and TF synthesis *rate* (molecules per hour) after to differentiation.
* maxesPre[nElem], maxesPost[nElem]: Maximal molecule counts for each TF (mRNA if `isProt` == 0) before and after differentiation
* meansPre[nElem], meansPost[nElem]: Average molecule counts for each TF (mRNA if `isProt` == 0) before and after differentiation
* ttarray[nElem]: Array of "differentiation times" for each species, using DHT of each individual TF (mRNA) trajectory
* degrees[nElem]: array of mRNA indegrees from `*.prot.degrees.txt`
