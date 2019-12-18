# Generalized Interaction Trees

Code for paper "Tree Based Subgroup Identification for Observational Data".

* `main/`: code to run simulations in the main manuscript.
  + `MainSimulationsCt.R`: run causal tree algorithms.
  + `MainSimulationsIpw.R`: run inverse probability weighting interaction trees.
  + `MainSimulationsG.R`: run G-computation interaction trees.
  + `MainSimulationsDr.R`: run doubly robust interaction trees.
  
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
  
  
  