# GRN Simulation (`mainsim`)

Note: This document and source are currently intended for reviewers of the associated miRNA-control manuscript/submission by R. Posner and R. Laubenbacher. This is not a "general-purpose" readme at the moment, but will be expanded accordingly following review and/or acceptance.

## Invocation

The simulator should be called as `mainsim -f grn_dir -n activeMirs -s nSimulations` where:
* `grn_dir` is a path to the directory containing the network in a file called "genes.json". For example, if the network generator is used, the path may look like `output/net.md.[MD5]/ngs.[networkGenSeed]`.
* `activeMirs` is the number of miRNAs from the network that should be activated. The specific activated miRNAs are randomly chosen for each invocation, but can be typically found in the final output directory (see below).
* `nSimulations` is the number of times to run the actual simulation

Note: a warning will be generated when `activeMirs` is invalid - the warning will indicate what was *actually* run

## Output

The simulator will create the following output:

1. Within `grn_dir`, a subdirectory called `miRs.[activeMirs]` will be created if it does not exist.
1. A further subdirectory will be created, called `states.[AABABABBAA]`, with 40 total A's/B's. This corresponds to the first 40 miRNA active states, where an "A" ("B") in position **x** indicates that the x-th miRNA is active (inactive, respectively). Although these strings are not unique (when more than 40 miRNAs are present), the odds of different active/inactive miRNA sequence producing the same directory name are about one per trillion.
1. Within the `states.[AB...]` subdirectory, if it does not already exist, a copy of `genes.json` will be placed within. This file has the inactive elements (and their connections) removed, and should be the same network for all output in its directory.
1. (Typically) three tab-delimited files will be created. They are:
  * `sim.[simSeed].prot.txt`
  * `sim.[simSeed].rna.txt` - These files contain the molecule/transcript counts for each gene at approximately uniform time intervals. Each row is \# TFs + 1 columns wide. The first column contains the time value (`clockTime`) at which the values were sampled. The remaining columns list copy numbers. There is a header row containing the names for each logged entity.
  * `sim.[simSeed].prot.degrees.txt` - The number of miRNAs targeting each (TF-encoding) mRNA species. Contains two rows. The first row has tab-delimited protein names, the second lists the (mRNA in-)degree.
Note: the number of files and the specific output can be modified using different compile-time-defined values - listed above is the "default" configuration.
Note 2: If one wants to regenerate the sim output, the shortened `genes.json` file placed in the simulation output directory and the value of `simSeed` *completely define* the output (a little cumbersome to recreate currently, as the `FIXSEED` macro can be used to fix the simulation seed, which is currently a compile-time constant).


## Quirks

The simulator has been optimized for performance on the computer used to perform the simulations (Intel i7-8700K model CPU), and some architecture-specific elements have been used. Some minor changes to the source have been made to improve multi-platform support. Specifically, in `include/globals.h`, the occupancy vector is represented using a 64-bit integer, and the macro `MAX_NCONNECTIONS` is accordingly set to 64. For the purposes of the paper, the Intel-specific  `_m256` intrinsic was used with `MAX_NCONNECTIONS` set to 255. Also, the copy numbers for coding and noncoding genes were hard-coded in the macros `CODING_QTY` and `NONCODING_QTY`. The maximum number of transcripts for any gene is defined using the macro `MAX_MRNA_QTY` and when compiled using `make opt`, it _performs no bounds checking_, so caution must be used when changing mRNA transcription or decay rates to ensure that this number is sufficiently high. Bounds checking takes place for other compilation options. These quirks are noted in the source code as well
