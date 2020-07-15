# Causal Interaction Trees

Code for paper "Causal Interaction Trees: Finding Subgroups with Heterogeneous Treatment Effects in Observational Data" ([arXiv link](https://arxiv.org/abs/2003.03042)).

* `main/`: code for simulations and data analysis in the main manuscript.
  * Simulations:
    + `CausalTreeOriginal.R`: run Original Causal Tree algorithms.
    + `CausalTreeBest.R`: run Best Causal Tree algorithms.
    + `MainIpw.R`: run Inverse Probability Weighting Causal Interaction Tree algorithms.
    + `MainG.R`: run G-formula Causal Interaction Tree algorithms.
    + `MainDr.R`: run Doubly Robust Causal Interaction Tree algorithms.
  * Analysis of the SUPPORT study:
    + `RhcOptCt.R`: run Causal Tree algorithms.
    + `RhcOptIpw.R`: run Inverse Probability Weighting Causal Interaction Tree algorithms.
    + `RhcOptG.R`: run G-formula Causal Interaction Tree algorithms.
    + `RhcOptDr.R`: run Doubly Robust Causal Interaction Tree algorithms.
    + `RhcRootBi.R`: calculate the bootstrap interval for the root node.
    + `RhcFrstSpltBi.R`: calculate the bootstrap intervals for the two subgroups given the first split in the large tree built by DR-CIT.
  * `Run*.R`: shell scripts to run associated R scripts on `slurm`.
  
* `Functions/`: functions to implement the Inverse Probability Weighting, G-computation, and Doubly Robust Interaction trees.
  + `CvMethod1.R`: functions for the main final tree selection with different estimators. 
  + `CvMethod2.R`: functions for the alternative final tree selection with different estimators. 
  + `EstCondEff.R`: `est.cond.eff()` fits the outcome model. `est.prop.sc()` fits the propensity score model.
  + `EstDrTempFunc.R`: splitting and evaluation functions for the doubly robust estimator.
  + `EstGTempFunc.R`: splitting and evaluation functions for the G-computation estimator.
  + `EstIpwTempFunc.R`: splitting and evaluation functions for the inverse probability weighting estimator.
  + `EvalMeas.R`: functions to evaluate final trees with the evaluation criteria described in paper.
  + `iTemp.R`: initialization function.
  + `library.R`: install required libraries, if not yet installed, and load the libraries.
  + `pruning.R`: build a maximum sized tree and create a sequence of candidate trees for final tree selection.
  + `SimData.R`: generate data under different settings.

* `seed1000.rda`: seeds for generating the data in simulations. The first 1000 seeds were used to generate data for the paper.

* `Appendix/`: Additional code to run simulations in the appendices.

* `Data/`: data for simulation results.
  
* `Results/`: code to generate figures and tables in the paper.
  + `main.R`: Figure 1 and Table 1 in the main manuscript.
  
  
  