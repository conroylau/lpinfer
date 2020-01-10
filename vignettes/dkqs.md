-   [Introduction](#introduction)
-   [Scope of the Vignette](#scope-of-the-vignette)
-   [Installation and Requirements](#installation-and-requirements)
-   [Usage Demonstration 1: Missing Data
    Problem](#usage-demonstration-1-missing-data-problem)
    -   [Background and Data](#data)
    -   [Syntax](#syntax)
    -   [Specifying the Functions](#functions)
    -   [Specifying the parameters](#specifying-the-parameters)
    -   [Output](#output)
-   [Help, Feature Requests and Bug
    Reports](#help-feature-requests-and-bug-reports)
-   [References](#references)

Introduction
------------

This package conducts inference on econometrics problems that can be
represented by quadratic programming. This package uses the
cone-tightening procedure of Deb et al. (2018) and the procedure of
choosing the optimal value of tau by Kamat (2019).

Scope of the Vignette
---------------------

This vignette is intended as a guide to using the `linearprog` package.
Readers may refer to section 4.2 of Deb et al. (2018) for details about
the cone-tightening procedure and the supplemental appendix of Kamat
(2019).

Installation and Requirements
-----------------------------

`linearprog` can be installed from our GitHub repository via

    devtools::install_github("conroylau/linearprog")

To use `linearprog`, one of the following packages for solving linear
and quadratic programs are required. There are four options for the
packages:

1.  Gurobi and the R package `gurobi` — Gurobi can be downloaded from
    [Gurobi Optimization](https://www.gurobi.com/). A Gurobi software
    license is required for this software, which can be obtained at no
    cost for academic researchers. The instructions for installing
    `gurobi` on R can be found
    [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html#r-package-installation).

2.  `limSolve` — a free and open-source package available on CRAN. This
    can be installed directly via the `install.packages` command in R.

3.  IBM ILOG CPLEX Optimization Studio (CPLEX) and one of the R packages
    below — CPLEX can be downloaded from
    [IBM](https://www.ibm.com/analytics/cplex-optimizer). The R packages
    are available on CRAN. Note that users who are using OS X would need
    to install the packages on the command line to link the package to
    the correct CPLEX library. Both R packages are open-source packages.
    Any one of the following packages can be used:

    1.  `Rcplex` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/Rcplex/INSTALL).

    2.  `cplexAPI` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/cplexAPI/INSTALL).

If no package is specified, one of the above packages will be
automatically chosen from those that are available.

Usage Demonstration 1: Missing Data Problem
-------------------------------------------

### Background and Data

This missing data problem is a classical problem due to Manski (1989),
which is demonstrated by the dataset below. The dataset contains 1,000
rows with 2 columns where `D` is binary and `Y` takes 11 distinct
values.

The sample data is included with the `linearprog` package. The sample
data has two columns and the following structure.

    library(linearprog)
    knitr::kable(head(sampledata, n = 10))

<table>
<thead>
<tr>
<th style="text-align:right;">
D
</th>
<th style="text-align:right;">
Y
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.9
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.9
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.0
</td>
</tr>
<tr>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.2
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.4
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.8
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.5
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.6
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.9
</td>
</tr>
<tr>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.7
</td>
</tr>
</tbody>
</table>

The two variables that are used are as follows:

-   `Y` is a multivariate discrete outcome variable taht takes value
    from 0 to 1 with step size 0.1.
-   `D` is a binary treatment where *D*<sub>*i*</sub> = 1 means that
    *Y*<sub>*i*</sub> is observed.

### Syntax

The main command for applying the cone-tightening procedure in
`linearprog` is called `dkqs_cone`.

    library(linearprog)
    dkqs_cone(df = sampledata, 
              A_obs = A_obs_twom, 
              A_tgt = A_target, 
              func_obs = func_two_moment, 
              beta_tgt = 0.375, 
              bs_seed = 1,
              bs_num = 100,
              p_sig = 2,
              tau_input = tau,
              solver = gurobi)

where

-   `df` refers to the data being used in the inference.
-   `A_obs` refers to the “observed matrix” in the inference.
-   `A_tgt` refers to the “target matrix” in the inference.
-   `func_obs` refers to the function that generates the vector of
    observed beta.
-   `beta_tgt` refers to the value of beta to be tested.
-   `bs_seed` refers to the starting value of the seed in bootstrap.
-   `bs_num` refers to the total number of bootstraps to be conducted.
-   `p_sig` refers to the number of decimal places in the p-value.
-   `tau_input` refers to the value of tau chosen by the user.
-   `solver` refers to the solver used to solve the linear and quadratic
    programs.

The following two sections explain how the required parameters and
functions can be constructed in order to apply the `dkqs_cone` function.

### Specifying the Functions

The `dkqs_cone` function provides flexibility for users to specify the
functions to generate the observed value of beta from the dataset. The
requirement for the function provided by the user is as follows:

-   It only takes one argument, which is the dataset.
-   It only returns a numeric vector that corresponds to the vector for
    the observed beta.

#### Full Information Approach

The following is an example of defining the function for the full
information approach:

    func_full_info <- function(df){
      beta = NULL
      y_list = sort(unique(df[,"Y"]))
      n = dim(df)[1]
      yn = length(y_list)
      for (i in 1:yn){
        beta_i = sum((df[,"Y"] == y_list[i]) * (df[,"D"] == 1))/n
        beta = c(beta,c(beta_i))
      }
      beta = as.matrix(beta)
      return(beta)
    }

The `func_full_info` function returns a vector of length *J* + 1 where
element *i* of the vector refers to the probability that the
corresponding value of *y*<sub>*i*</sub> is observed.

#### Two Moments Approach

The following is an example of defining the function for the two moments
approach:

    func_two_moment <- function(df){
      beta = matrix(c(0,0), nrow = 2)
      n = dim(df)[1]
      beta[1] = sum(df[,"Y"] * df[,"D"])/n
      beta[2] = sum(df[,"D"])/n
      return(beta)
    }

The `func_two_moment` function returns a vector with two elements that
corresponds to the two moments.

### Specifying the parameters

As shown in the syntax section, the matrices `A_obs` and `A_tgt` have to
be defined in order to use the function. To construct the two matrices,
the following parameters are needed:

    N = dim(sampledata)[1]
    J1 = length(unique(sampledata[,"Y"]))
    yp = seq(0,1,1/(J1-1))

With the above quantities, the two matrices can be defined as follows:

    A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                    byrow = TRUE)
    A_target = matrix(c(yp, yp), nrow = 1)

Lastly, the value of tau can be defined freely by the user as long as
the quadratic program is feasible. Here, we choose the value of tau
based on the formula from the supplemental appendix in page 15 of Kamat
(2019):

    tau = sqrt(log(N)/N)

### Output

The following are the output when the two moments approach is used with
the `gurobi` solver to test the hypothesis that `beta_tgt` is 0.375.

    library(linearprog)
    dkqs_cone(df = sampledata, 
             A_obs = A_obs_twom, 
             A_tgt = A_target, 
             func_obs = func_two_moment, 
             beta_tgt = 0.375, 
             bs_seed = 1,
             bs_num = 100,
             p_sig = 3,
             tau_input = tau,
             solver = "gurobi")
    #> Linear and quadratic programming solver used: gurobi.
    #> ----------------------------------- 
    #> Test statistic: 0.06724.
    #> p-value: 0.253.
    #> Value of tau used: 0.08311.

Help, Feature Requests and Bug Reports
--------------------------------------

Please post an issue on the [GitHub
repository](https://github.com/conroylau/linearprog/issues).

References
----------

Deb, R., Y. Kitamura, J. K. H. Quah, and Stoye J. 2018. “Revealed Price
Preference: Theory and Empirical Analysis.” *Working Paper*.

Kamat, V. 2019. “Identification with Latent Choice Sets.” *Working
Paper*.

Manski, C. F. 1989. “Anatomy of the Selection Problem.” *The Journal of
Human Resources* 24: 343–60.
