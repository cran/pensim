  <!-- badges: start -->
  [![R-CMD-check](https://github.com/waldronlab/pensim/workflows/R-CMD-check/badge.svg)](https://github.com/waldronlab/pensim/actions)
  [![Coverage Status](https://codecov.io/github/waldronlab/pensim/coverage.svg?branch=master)](https://codecov.io/github/waldronlab/pensim?branch=master)
  [![](https://cranlogs.r-pkg.org/badges/pensim)](https://cran.r-project.org/package=pensim)
  <!-- badges: end -->

# pensim: Simulation of High-dimensional Data and Parallelized Repeated Penalized Regression

Waldron L, Pintilie M, Tsao M-S, Shepherd FA, Huttenhower C, Jurisica
I: [Optimized application of penalized regression methods to diverse
genomic data.][manuscript] *Bioinformatics* 2011,
27:3399–3406. 

# Abstract

**Motivation**: Penalized regression methods have been adopted widely
  for high-dimensional feature selection and prediction in many
  bioinformatic and biostatistical contexts. While their theoretical
  properties are well-understood, specific methodology for their
  optimal application to genomic data has not been determined.

**Results**: Through simulation of contrasting scenarios of correlated
  high-dimensional survival data, we compared the LASSO, Ridge and
  Elastic Net penalties for prediction and variable selection. We
  found that a 2D tuning of the Elastic Net penalties was necessary to
  avoid mimicking the performance of LASSO or Ridge
  regression. Furthermore, we found that in a simulated scenario
  favoring the LASSO penalty, a univariate pre-filter made the Elastic
  Net behave more like Ridge regression, which was detrimental to
  prediction performance. We demonstrate the real-life application of
  these methods to predicting the survival of cancer patients from
  microarray data, and to classification of obese and lean individuals
  from metagenomic data. Based on these results, we provide an
  optimized set of guidelines for the application of penalized
  regression for reproducible class comparison and prediction with
  genomic data.

**Availability and Implementation**: A parallelized implementation of
  the methods presented for regression and for simulation of synthetic
  data is provided as the pensim R package, available at
  http://cran.r-project.org/web/packages/pensim/index.html.

# Figure 1

<a href="https://bioinformatics.oxfordjournals.org/content/27/24/3399"><img src="https://github.com/waldronlab/schematics/raw/master/jpgs/F1_pensim.jpg"/></a>


**(A) Methodology for model selection and validation of high-dimensional
data.** Objectives include both feature selection and outcome
prediction, e.g. for patient survival given tumor gene expression
data. A nearly unbiased assessment of prediction accuracy for small
samples sizes is obtained by repeating all steps of model selection in
each iteration of the cross-validation. Variable selection and model
conditioning are achieved within the training sets by an optional,
permissive univariate pre-filter followed by repeated cross-validation
for parameter tuning. These steps are detailed in [Section 4]. **(B)
Over-fitting occurs in spite of tuning the models by cross-validation,
as evidenced by reduced prediction accuracy in simulated test sets
compared to resubstitution of training data.**

[Section 4]: https://academic.oup.com/bioinformatics/article/27/24/3399/306905#SEC4
[manuscript]: https://doi.org/10.1093/bioinformatics/btr591
