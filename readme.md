# EcoPE

This repo contains code and data to support Terry (2026) *Hierarchical permutation entropy regression identifies decadal decline in the intrinsic predictability of butterfly populations*

You may reuse the *new* code and data as you see fit (per MIT License.) 

For reproducibility only, original open datasets from the UKBMS, OS and Met Office are bundled raw in this repository. Any further reuse of that data must refer to the original sources and follow those original (open) licensing and acknowledgement notices. 

## Contents:

Core code and results are in the numbered R markdown files as bth raw .rmd and knit .htmls. See session information at the end of each html for exact package versions used. 


- 1-3 Explore increasingly complex regression models for permutation entropy in simulated systems
- 4 Filters and processes the original UKBMS butterfly data into a motif format. 
- 5a Fits the PE model to the butterfly data
- 5b Explores the model fits and conducts cross-species analyses
- 6 Conducts a parallel analysis of the butterfly data using auto-regressive models. 


Substantial model outputs are provided in `model_fits/` and `Butterfly_fits/`so that intensive models do not need to be refitted.


Stan and R function scripts are in their respective folders. Stan models are heavily annotated and have been developed for **`CmdStan` 2.39** - NB certain features will not work successfully with earlier versions of Stan, and standard `RStan` releases can lag. 









