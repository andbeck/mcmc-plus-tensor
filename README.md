mcmc-plus-tensor
================

R code, simulated data etc for Robinson and Beckerman 2013 Ecology Letters doi: 10.1111/ele.12047.

Bayesian MCMC based Matrix Comparison and Tensor Based Analysis of Phenotypic Plasticity, Ageing and Population Gradients too.

Updated 11 May 2016
1) unlocked restriction on 1000 samples.  Tools were requiring >=1000, and using 1000 as sample size.  Now uses sample size reported for model from MCMCglmm, which is the dimension of the $VCV output.

2) removed deprecated EISPACK call in eigen usage

3) Cleaned up reporting to map onto new sample size