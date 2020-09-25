lpinfer: An R Package for Inference in Linear Programs
================
Conroy Lau and Alexander Torgovitsky

Introduction
------------

This package provides a set of methods for estimation and statistical inference on the solutions (optimal values) of linear programs.
The motivation is partially identified models.
Identified sets in many partially identified models can be characterized by the solutions to two linear programs (one min, one max), with the interval between comprising the whole identified set.
The module is designed with this situation as the primary use-case, although it might be of interest in other stochastic programming applications where nonparametric inference is desired.

Installation
-----------------------------

`lpinfer` can be installed from its GitHub repository via

```r
devtools::install_github("conroylau/lpinfer")
```

Using `lpinfer` requires having a package for solving linear and quadratic programs.
The following options are supported:

1. **(Strongly recommended)** Gurobi and the R package `gurobi` — Gurobi can be
   downloaded from
    [Gurobi Optimization](https://www.gurobi.com/). A Gurobi software
    license is required. The license can be obtained at no cost for
    academic researchers. [This guide](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html#r-package-installation)
    provides a very clear set of instructions for installing `gurobi` in R.

2. IBM ILOG CPLEX Optimization Studio (CPLEX) and one of the R packages
    below — CPLEX can be downloaded from
    [IBM](https://www.ibm.com/analytics/cplex-optimizer). A CPLEX
    software license is required, which can be obtained at no cost for
    academic researchers. There are two free and open-source R packages
    that provide APIs to CPLEX, and `lpinfer` supports both:

    1. `Rcplex` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/Rcplex/INSTALL).

    2. `cplexAPI` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/cplexAPI/INSTALL).

3. `limSolve`, a free and open-source package available on CRAN.

4. `lpSolveAPI`, another free and open-source package available on CRAN.
This package cannot solve quadratic programs, so cannot be used with some of the methods in `lpinfer`.

Current Features
-----------------------------

- Estimation using the procedure of [Mogstad, Santos and Torgovitsky (2018, _Econometrica_)](https://doi.org/10.3982/ECTA15463).
- Statistical inference (testing and confidence intervals) using:
    - The obtuse angle procedure developed by [Fang, Santos, Shaikh and Torgovitsky (2020, working paper)](https://a-torgovitsky.github.io/fsst.pdf).
    - The cone-tightening procedure developed by [Deb, Quah, Kitamura and Stoye (2018, working paper)](https://arxiv.org/abs/1801.02702v2).
    - Profiled subsampling, as proposed by [Romano and Shaikh (2008, _Journal of Statistical Planning and Inference_)](https://doi.org/10.1016/j.jspi.2008.03.015).
    - The direct bootstrap procedure proposed by [Cho and Russell (2019, working paper)](https://arxiv.org/abs/1810.03180).
- Support for parallel programming using the `future` and `futute.apply` packages.
- **(experimental)** Automatic conversion from extensive-form linear programs to standard form.

Planned Features
-----------------------------

More inference procedures will be added soon.

Usage
-----------------------------
A detailed vignette is available [here](./vignette/vignette.pdf).

Help, Feature Requests, Bug Reports, and Contributing
--------------------------------------

Please post an issue on the [issues
page](https://github.com/conroylau/lpinfer/issues) of the GitHub
repository.

If you have developed a procedure that you would like implemented, please
contact us. We may be able to work together to include your method.

How to Cite
----------
If you use `lpinfer` in your research please cite the Zenodo DOI.
