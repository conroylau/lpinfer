linearprog: An R Package for Inference in Linear and Quadratic Programs
================
-   [Introduction](#introduction)
-   [Scope of the Vignette](#scope-of-the-vignette)
-   [Installation and Requirements](#installation-and-requirements)
-   [Usage Demonstration 1: Missing Data
    Problem](#usage-demonstration-1-missing-data-problem)
    -   [Background and Data](#data)
    -   [Specifying the Functions](#functions)
    -   [Specifying the Parameters](#specifying-the-parameters)
    -   [Calculating *p*-value](#calculating-p-value)
    -   [Constructing Confidence
        Intervals](#constructing-confidence-intervals)
    -   [Constructing Multiple Confidence
        Intervals](#constructing-multiple-confidence-intervals)
-   [Help, Feature Requests and Bug
    Reports](#help-feature-requests-and-bug-reports)
-   [References](#references)

Introduction
------------

This package conducts inference on econometrics problems that can be
studied by linear and quadratic programming using the cone-tightening
procedure of Deb et al. (2018).

Scope of the Vignette
---------------------

This vignette is intended as a guide to use the `linearprog` package.
Readers may refer to section 4.2 of Deb et al. (2018) for details about
the cone-tightening procedure and the supplemental appendix of Kamat
(2019).

Installation and Requirements
-----------------------------

`linearprog` can be installed from our GitHub repository via

    devtools::install_github("conroylau/linearprog")

To use `linearprog`, one of the following packages for solving linear
and quadratic programs is required. There are four options for the
solver:

1.  Gurobi and the R package `gurobi` — Gurobi can be downloaded from
    [Gurobi Optimization](https://www.gurobi.com/). A Gurobi software
    license is required, which can be obtained at no cost for academic
    researchers. The instructions for installing `gurobi` on R can be
    found
    [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html#r-package-installation).

2.  IBM ILOG CPLEX Optimization Studio (CPLEX) and one of the R packages
    below — CPLEX can be downloaded from
    [IBM](https://www.ibm.com/analytics/cplex-optimizer). A CPLEX
    software license is required, which can be obtained at no cost for
    academic researchers. There are two open-source and free R packages
    that uses CPLEX, and users are free to choose one of them. In
    addition, both packages have to be installed on the command line to
    link the package to the correct CPLEX library. The two packages’
    name and installation instrutions are as follows:

    1.  `Rcplex` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/Rcplex/INSTALL).

    2.  `cplexAPI` — the instructions to install the R package can be
        found
        [here](https://cran.r-project.org/web/packages/cplexAPI/INSTALL).

3.  `limSolve` — a free and open-source package available on CRAN. This
    can be installed directly via the `install.packages` command in R.

If no package is specified, one of the above packages will be
automatically chosen from those that are available.

Usage Demonstration 1: Missing Data Problem
-------------------------------------------

### Background and Data

One application of this `linearprog` package is to study the classical
missing data problem due to Manski (1989), where the sharp bounds for
the identified set of the expected value of the observed value can be
constructed by linear and quadratic programming. The `linearprog`
package can be used to construct the *p*-value and confidence interval
for the missing data problem.

As an illustration, the dataset below studies the missing data problem
and contains 1,000 simulated data with 2 columns. This dataset is
included in the `linearprog` package.
``` r
library(linearprog)
knitr::kable(head(sampledata, n = 10))
``` 
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

where

-   `Y` is a multivariate discrete outcome variable that takes value
    from 0 to 1 with step size 0.1.
-   `D` is a binary treatment where *Y*<sub>*i*</sub> is observed for
    *D*<sub>*i*</sub> = 1 and not observed for *D*<sub>*i*</sub> = 0.

### Specifying the Functions

To conduct the tests in this `linearprog` package, users need to specify
the functions to generate the observed value of beta from the dataset.
This package provides the flexibility for users to specify their own
function to generate the value of observed value of beta. The
requirement for the function provided by the user is as follows:

-   The function’s only argument is the dataset.
-   The function only returns a numeric vector that corresponds to the
    vector for the observed beta.

#### Full Information Approach

The following is an example of defining the function for the full
information approach:
``` r
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
```
The `func_full_info` function returns a vector of length that is equal
to the number of distinct observations for *Y* where element *i* of the
vector refers to the probability that the corresponding value of
*y*<sub>*i*</sub> is observed.

#### Two Moments Approach

The following is an example of defining the function for the two moments
approach:
``` r
func_two_moment <- function(df){
  beta = matrix(c(0,0), nrow = 2)
  n = dim(df)[1]
  beta[1] = sum(df[,"Y"] * df[,"D"])/n
  beta[2] = sum(df[,"D"])/n
  return(beta)
}
```
The `func_two_moment` function returns a vector with two elements that
corresponds to the two moments **E**[*Y*<sub>*i*</sub>] and
**E**[*Y*<sub>*i*</sub>*D*<sub>*i*</sub>].

### Specifying the Parameters

To conduct the inference in this package, the matrices `A_obs` and
`A_tgt` have to be defined in order to use the function. To construct
the two matrices, the following parameters are needed:
``` r
N = dim(sampledata)[1]
J1 = length(unique(sampledata[,"Y"]))
yp = seq(0,1,1/(J1-1))
```
With the above quantities, the following matrices can be defined as
follows:
``` r
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                byrow = TRUE)
A_target = matrix(c(yp, yp), nrow = 1)
```
The matrix `A_obs_twom` refers to the observed matrix for the two
moments approach. If users would prefer using the full information
approach, the following matrix that correspond to the full information
approach has to be defined:
``` r
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
```
Lastly, the value of tau can be defined freely by the user as long as
the quadratic program is feasible. Here, we choose the value of tau
based on the formula from page 15 of the supplemental appendix of Kamat
(2019):
``` r
tau = sqrt(log(N)/N)
```
### Calculating *p*-value

The command for applying the cone-tightening procedure in `linearprog`
is called `dkqs_cone`. It takes the data, parameters and the function
that can calculate the observed value of beta and returns the *p*-value.

#### Syntax

This `dkqs_cone` command has the following syntax:
``` r
dkqs_cone(df = sampledata, 
          A_obs = A_obs_twom, 
          A_tgt = A_target, 
          func_obs = func_two_moment, 
          beta_tgt = 0.375, 
          bs_seed = 1,
          bs_num = 100,
          p_sig = 2,
          tau_input = tau,
          solver = "gurobi",
          cores = 1,
          progress = TRUE)
```
where

-   `df` refers to the data set.
-   `A_obs` refers to the “observed matrix”.
-   `A_tgt` refers to the “target matrix”.
-   `func_obs` refers to the function that generates the vector of
    observed beta.
-   `beta_tgt` refers to the value of beta to be tested.
-   `bs_seed` refers to the starting value of the seed in bootstrap.
-   `bs_num` refers to the total number of bootstraps to be conducted.
-   `p_sig` refers to the number of decimal places in the *p*-value.
-   `tau_input` refers to the value of tau chosen by the user.
-   `solver` refers to the name of the solver used to solve the linear
    and quadratic programs.
-   `cores` refers to the number of cores to be used in the parallelized
    for-loop for computing the bootstrap test statistics. See <span
    id="parallel-dkqs">here</span> for more details.
-   `progress` refers to the boolean variable for whether the result
    messages should be displayed in the procedure of calculating the
    *p*-value.

#### Using Parallel Programming

The `linearprog` package supports the use of parallel programming in
computing the bootstrap test statistics to reduce the computational
time. To use parallel programming, specify the number of cores that you
would like to use in the argument `cores`. If you do not want to use
parallel programming, you may input 1 or any other non-numeric variables
in `cores`.

For best performance, it is advisible to specify the number of cores to
be less than or equal to the cores that you have on your machine.

#### Output

The followings are the output when the **two moments approach** is used
with the `gurobi` solver to test the hypothesis that `beta_tgt` is
0.375.
``` r
dkqs_cone(df = sampledata, 
         A_obs = A_obs_twom, 
         A_tgt = A_target, 
         func_obs = func_two_moment, 
         beta_tgt = 0.375, 
         bs_seed = 1,
         bs_num = 100,
         p_sig = 3,
         tau_input = tau,
         solver = "gurobi",
         cores = 8,
         progress = TRUE)
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used in bootstrap procedure: 8.
#> ----------------------------------- 
#> Test statistic: 0.06724.
#> p-value: 0.253.
#> Value of tau used: 0.08311.
``` 
Alternatively, the followings are the output when the **full information
approach** is used with the `gurobi` solver to test the hypothesis that
`beta_tgt` is 0.375.
``` r
dkqs_p = dkqs_cone(df = sampledata, 
                   A_obs = A_obs_full, 
                   A_tgt = A_target, 
                   func_obs = func_full_info, 
                   beta_tgt = 0.375, 
                   bs_seed = 1,
                   bs_num = 100,
                   p_sig = 3,
                   tau_input = tau,
                   solver = "gurobi",
                   cores = 8,
                   progress = TRUE)
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used in bootstrap procedure: 8.
#> ----------------------------------- 
#> Test statistic: 0.01746.
#> p-value: 0.253.
#> Value of tau used: 0.08311.
```
The results from the two approach should give the same *p*-value while
the test statistic will be different as different moments are considered
in the problem.

The message generated at the end of the function can be re-generated by
`print(dkqs_p)`. All message displayed when the function is calculating
the *p*-value can be re-generated by `summary(dkqs_p)`.

### Constructing Confidence Intervals

The command for constructing the two-sided confidence interval for a
statistical test in the `linearprog` packge is called `qpci`. The
confidence interval is constructed by evaluating the *p*-value of a test
and applying the biscetion method.

#### Syntax

The syntax of the `qpci` function is as follows:
``` r
qpci(f = dkqs_cone, 
     farg = dkqs_farg, 
     alpha = 0.05, 
     lb0 = NULL, 
     lb1 = NULL, 
     ub0 = NULL, 
     ub1 = NULL, 
     tol = 0.0001, 
     max_iter = 20, 
     df_ci = NULL,
     progress = FALSE)
``` 
where

-   `f` refers to the function that represents a testing procedure.
-   `farg` refers to the list of arguments to be passed to the function
    of testing procedure.
-   `alpha` refers to the signifance level of the test.
-   `lb0` refers to the logical lower bound for the confidence interval.
-   `lb1` refers to the maximum possible lower bound for the confidence
    interval.
-   `ub0` refers to the logical upper bound for the confidence interval.
-   `ub1` refers to the minimum possible upper bound for the confidence
    interval.
-   `tol` refers to the tolerance level in the bisection method.
-   `max_iter` refers to the maximum number of iterations in the
    bisection method.
-   `df_ci` refers to dataframe that consists of the points and the
    corresponding *p*-values that have been tested in constructing the
    confidence intervals.
-   `progress` refers to the boolean variable for whether the result
    messages should be displayed in the procedure of constructing
    confidence interval.

#### Specifying the Argument

To use the `qpci` function, the arguments for the function of the test
statistic has to be specified and passed to `farg`. For instance, if the
test `dkqs_cone` is used, then the arguments can be defined as follows:
``` r
dkqs_farg = list(df = sampledata,
                 A_obs = A_obs_full,
                 A_tgt = A_target,
                 func_obs = func_full_info,
                 bs_seed = 1,
                 bs_num = 100,
                 p_sig = 2,
                 tau_input = tau,
                 solver = "gurobi",
                 cores = 8,
                 progress = FALSE)
``` 
Note that the argument for the target value of beta, i.e. the value to
be tested under the null, is not required in the above argument
assignment.

#### Specifying the Data Frame `df_ci`

If the *p*-values at certain points have already been evaluated, users
can store them in a data frame and pass it to the function `qpci`. The
requirement for the data frame is as follows:

-   The data frame can only has two columns. The first column is `point`
    (which contains the values of betas that has been evaluated) and the
    second column is `value` (which corresponds to the *p*-values being
    evaluated).
-   The data frame can only contain numeric values.

#### Output

The following shows a sample output of the function `qpci` that is used
to the confidence interval for the test `dkqs_cone` with significance
level 0.05.
``` r
qpci_dkqs = qpci(f = dkqs_cone, 
                 farg = dkqs_farg, 
                 alpha = 0.05, 
                 lb0 = 0, 
                 lb1 = 0.4, 
                 ub0 = 1, 
                 ub1 = 0.6, 
                 tol = 0.001, 
                 max_iter = 5, 
                 df_ci = NULL, 
                 progress = TRUE)
#> 
#> === Computing upper bound of confidence interval ===
#> >>> Evaluating the first left end-point
#>       * Point being evaluated: 0.6
#>       * p-value: 0.78
#>       * Decision: Do not reject
#> >>> Evaluating the first right end-point
#>       * Point being evaluated: 1
#>       * p-value: 0
#>       * Decision: Reject
#> >>> Iteration 1
#>       * Point being evaluated: 0.8
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.6, 1]
#> >>> Iteration 2
#>       * Point being evaluated: 0.7
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.6, 0.8]
#> >>> Iteration 3
#>       * Point being evaluated: 0.65
#>       * p-value: 0.05
#>       * Decision: Do not reject
#>       * Current interval: [0.6, 0.7]
#> >>> Iteration 4
#>       * Point being evaluated: 0.675
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.65, 0.7]
#> >>> Iteration 5
#>       * Point being evaluated: 0.662
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.65, 0.675]
#>       * Reached the maximum number of iterations.
#> 
#> === Computing lower bound of confidence interval ===
#> >>> Evaluating the first left end-point
#>       * Point being evaluated: 0
#>       * p-value: 0
#>       * Decision: Reject
#> >>> Evaluating the first right end-point
#>       * Point being evaluated: 0.4
#>       * p-value: 0.79
#>       * Decision: Do not reject
#> >>> Iteration 1
#>       * Point being evaluated: 0.2
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0, 0.4]
#> >>> Iteration 2
#>       * Point being evaluated: 0.3
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.2, 0.4]
#> >>> Iteration 3
#>       * Point being evaluated: 0.35
#>       * p-value: 0
#>       * Decision: Reject
#>       * Current interval: [0.3, 0.4]
#> >>> Iteration 4
#>       * Point being evaluated: 0.375
#>       * p-value: 0.26
#>       * Decision: Do not reject
#>       * Current interval: [0.35, 0.4]
#> >>> Iteration 5
#>       * Point being evaluated: 0.363
#>       * p-value: 0.03
#>       * Decision: Do not reject
#>       * Current interval: [0.35, 0.375]
#>       * Reached the maximum number of iterations.
#> -----------------------------------
#> Total number of iterations: 10.
#> Tolerance level: 0.001.
#> Confidence interval: [0.356, 0.656].
```
The message generated at the end of the function can be re-generated by
`print(qpci_dkqs)`. All message displayed when the function is
constructing the confidence interval can be re-generated by
`summary(qpci_dkqs)`.

### Constructing Multiple Confidence Intervals

The `many_qpci` function is a wrapper for the `qpci` function where the
user can pass multiple values of alpha so multiple confidence intervals
can be constructed in one command.

#### Syntax

The syntax for the `many_qpci` function is as follows:
``` r
many_qpci(f = dkqs_cone, 
          farg = dkqs_farg, 
          alphas = c(0.01, 0.05), 
          lb0 = NULL, 
          lb1 = NULL, 
          ub0 = NULL, 
          ub1 = NULL, 
          tol = 0.0001, 
          max_iter = 20, 
          df_ci = NULL,
          progress_one = TRUE,
          progress_many = TRUE)
```
where

-   `alphas` refers to the list of significance levels to be used in
    constructing the confidence intervals.
-   `progress_one` refers to the boolean variable for whether the result
    messages should be displayed in running the function `qpci`.
-   `progress_many` refers to the boolean variable for whether the
    result messages should be displayed in running the function
    `many_qpci`.
-   The rest of the arguments are the same as the function `qpci` and
    they can be found [here](#qpci_syntax).

#### Output

The following shows a sample output of the function `qpci_many` that is
used to the confidence intervals for the test `dkqs_cone` with
significance level 0.01, 0.05 and 0.1.
``` r
many_qpci(f = dkqs_cone, 
          farg = dkqs_farg, 
          alphas = c(0.01, 0.05, 0.1), 
          lb0 = 0, 
          lb1 = 0.4, 
          ub0 = 1, 
          ub1 = 0.6, 
          tol = 0.001, 
          max_iter = 10, 
          df_ci = NULL, 
          progress_one = FALSE, 
          progress_many = TRUE)
#> Confidence interval for significance level 0.01: [0.35742, 0.65977].
#> Confidence interval for significance level 0.05: [0.36211, 0.65352].
#> Confidence interval for significance level 0.1: [0.36523, 0.65039].
``` 
Help, Feature Requests and Bug Reports
--------------------------------------

Please post an issue on the [GitHub
repository](https://github.com/conroylau/linearprog/issues), and we are
happy to help.

References
----------

Deb, R., Y. Kitamura, J. K. H. Quah, and Stoye J. 2018. “Revealed Price
Preference: Theory and Empirical Analysis.” *Working Paper*.

Kamat, V. 2019. “Identification with Latent Choice Sets.” *Working
Paper*.

Manski, C. F. 1989. “Anatomy of the Selection Problem.” *The Journal of
Human Resources* 24: 343–60.
