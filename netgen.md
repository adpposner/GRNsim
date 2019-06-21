# Network generation

Note: This document and source are currently intended for reviewers of the associated miRNA-control manuscript/submission by R. Posner and R. Laubenbacher. This is not a "general-purpose" readme at the moment, but will be expanded accordingly following review and/or acceptance.


Gene regulatory networks (GRNs) are generated from a set of pre-specified parameters. The parameters are organized topically in the file *config.xml*, and are listed below.

## Invocation
`config.xml` must be in same directory as network generation program at time of invocation, and will be read automatically. Therefore, one only needs to call `mainnetgen` to produce networks. Additionally, the base directory must contain a subdirectory called `output` to write networks to.

## Output
`mainnetgen` will produce an MD5 hash of *all parameter values, except those in the `randomSeeds` block* from `config.xml`. Then:
1. If it does not already exist, the folder `output/net.md.[MD5]` will be created, and a copy of `config.xml` will be placed within.
    * Hence, parameters used for network generation will be stored with the networks themselves. The MD5 hash is used to identify networks created with *identically distributed* parameters.
1. For each invocation of `mainnetgen` with a particular value of `networkGenSeed`, a subdirectory of  `output/net.md.[MD5]` will be created called `ngs.[networkGenSeed]`. 
1. Within that "ngs" subdirectory, the generated network will be written to the file `genes.json.` The network is *completely* specified by the `config.xml` file in  its parent directory, and the value in the directory name (i.e. `ngs.[networkGenSeed]`).
See `GRN input description` for details.




## GRN parameters (config.xml)

* globalDimensions
  * nMess - \# of genes (DNA) encoding TFs/mRNAs/TF species
  * nMicro - \# of DNA elements encoding miRNAs/miRNA species
* connectionParameters
  * miRNAs
    * inOutDegree - (value = "in" or "out") specifies whether values correspond to in or out degrees for miRNA/mRNA connections (directions with respect to miRNAs)
    * degree_low - minimum degree for miRNA/mRNA connections
    * degree_high - maximum degree for miRNA/mRNA connections
    * effectProb - probability of a "connection" - the (in or out) degree distributions are given by *deg_low + Binomial(deg_high-deg_low, effectProb)*
  * TFs - Separated by type of DNA species with which TF interacts - experimental
    * coding
      * Same parameters/names as for (miRNAs) - correspond to TF/TF-encoding DNA interactions
    * noncoding (experimental, for manuscript only the "coding" is used, and it applies to all TF-DNA interactions, both coding and noncoding)
      * Same parameters/names as for (miRNAs) - correspond to TF/miRNA-encoding DNA interactions
* randomSeeds
  * networkGenSeed - Seed value given to PRNG for network generation. If all other parameters are same and networkGenSeed is same, identical outputs will be produced

* baseRates - reaction rates, same (except for "effect," for all genes/species)
  * mess
    * production - base mRNA transcription rate
    * decay - base mRNA decay rate
  * TF
    * production - base TF translation rate from *one* mRNA transcript molecule
    * decay - base TF decay rate
    * affinity - "tightness" of TF-DNA bond (see note below)
    * frequency - "likelihood" of TF-DNA binding (see note below)
    * effectMin - minimum effect of TF on transcription rate (<1 is inhibitory, >1 is activating)
    * effectMax - maximum effect of TF on transcription rate - effects are uniformly dist'd between effectMin and effectMax
  * micro
    * production - base miRNA production rate
    * decay - base miRNA decay rate
    * affinity - "tightness" of miRNA-mRNA bond (see note below)
    * frequency - "likelihood" of miRNA-mRNA binding (see note below)
    * effectStrength - multiplier effect of miRNA on translation - since we assume miRNAs are always inhibitory in manuscript/submission, we just use a fixed small constant (or 0 for complete abortion of translation)

Note: To obtain traditional "k_a" and "k_d" values, the formula k_a = 2\*A\*F / (1+A) and k_d = 2\*F / (1+A), where A is the "affinity" value, and F is the "frequency" value

  * defaultRate - Used only during development/debugging as a "catch-all" - typically use a negative value to cause errors when rates not properly specified elsewhere

* initialQuantities
  * coding - Copy \# for each TF-encoding gene
  * noncoding - Copy \# for each miRNA-encoding "gene"
  * mess - \# of initial mRNA transcripts for each gene
  * micro - \# of initial active miRNAs for each "gene"
  * TF - \# of initial TF molecules for each TF species

* prices (experimental, not used in manuscript/submission, but parsed by parameters.c and code for tracking/logging costs is already present and functional). It should be noted that the prices, as implemented, are currently constant and static, and do not account for important factors such as sequence length.
  * proteinCost - "energy" cost of producing one TF molecule from one mRNA transcript
  * mRNACost - "energy" cost of producing one mRNA transcript from DNA
  * miRNACost - "energy" cost of producing one mature miRNA molecule from DNA
