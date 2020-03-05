linearprog: An R Package for Inference in Linear and Quadratic Programs
================

-   [Introduction](#introduction)
-   [Scope of the Vignette](#scope-of-the-vignette)
-   [Installation and Requirements](#installation-and-requirements)
-   [Example Data](#data)
-   [Test 1: Cone-tightening
    procedure](#test-1-cone-tightening-procedure)
-   [Test 2: Subsampling](#test-2-subsampling)
-   [Constructing Confidence
    Intervals](#constructing-confidence-intervals)
    -   [Constructing a Single Confidence
        Interval](#constructing-a-single-confidence-interval)
    -   [Constructing Multiple Confidence
        Intervals](#constructing-multiple-confidence-intervals)
-   [Constructing Bounds subject to Shape
    Constraints](#constructing-bounds-subject-to-shape-constraints)
-   [Help, Feature Requests and Bug
    Reports](#help-feature-requests-and-bug-reports)
-   [References](#references)

Introduction
------------

This package conducts inference on econometrics problems that can be
studied by linear and quadratic programs. Currently, this package
supports the following tests:

1.  Cone-tightening procedure by Deb et al. (2018)
2.  Subsampling procedure

This package can compute the *p*-value based on the test, construct
confidence intervals and estimate the bounds subject to shape
constraints.

Scope of the Vignette
---------------------

This vignette is intended as a guide to use the `linearprog` package.
For the cone-tightening procedure, readers may refer to section 4.2 of
Deb et al. (2018) for details about the procedure and the supplemental
appendix of Kamat (2019).

Installation and Requirements
-----------------------------

`linearprog` can be installed from our GitHub repository via
```r
devtools::install_github("conroylau/linearprog")
```
To use most of the functions in `linearprog`, one of the following
packages for solving linear and quadratic programs is required. There
are four options for the solver:

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

The package `lpSolveAPI` is only supported in the function `estbounds`
when the L1-norm is used. This is a free and open-source package
available on CRAN. This can be installed directly via the
`install.packages` command in R.

Example Data
------------

The classical missing data problem due to Manski (1989) is used as an
example throughout the file to demomnstrate the functions in the
`linearprog` package. The problem is used because the sharp bounds for
the identified set of the expected value of the observed value can be
constructed by linear and quadratic programs.

As an illustration, the dataset below studies the missing data problem
and contains 1,000 simulated data with 2 columns. This dataset is
included in the `linearprog` package.
```r
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

Test 1: Cone-tightening procedure
---------------------------------

The function `dkqs` carries out the cone-tightening procedure that is
proposed by Deb et al. (2018).

#### Specifying the Functions

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
```r
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
```r
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

#### Specifying the Parameters

To conduct the inference in this package, the matrices `A_obs` and
`A_tgt` have to be defined in order to use the function. To construct
the two matrices, the following parameters are needed:
```r
N = dim(sampledata)[1]
J1 = length(unique(sampledata[,"Y"]))
yp = seq(0,1,1/(J1-1))
```
With the above quantities, the following matrices can be defined as
follows:
```r
A_obs_twom = matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
             byrow = TRUE)
A_target = matrix(c(yp, yp), nrow = 1)
```
The matrix `A_obs_twom` refers to the observed matrix for the two
moments approach. If users would prefer using the full information
approach, the following matrix that correspond to the full information
approach has to be defined:
```r
A_obs_full = cbind(matrix(rep(0,J1*J1), nrow = J1), diag(1, J1))
```
Lastly, the value of tau can be defined freely by the user as long as
the quadratic program is feasible. Here, we choose the value of tau
based on the formula from page 15 of the supplemental appendix of Kamat
(2019):
```r
tau = sqrt(log(N)/N)
```
#### Calculating *p*-value

The command for applying the cone-tightening procedure in `linearprog`
is called `dkqs`. It takes the data, parameters and the function that
can calculate the observed value of beta and returns the *p*-value.

#### Syntax

The `dkqs` command has the following syntax:
```r
dkqs(df = sampledata, 
     A_obs = A_obs_twom, 
     A_tgt = A_target, 
     func_obs = func_two_moment, 
     beta_tgt = 0.375, 
     bs_seed = 1,
     bs_num = 100,
     p_sig = 2,
     tau_input = tau,
     solver = gurobi,
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
```r
dkqs(df = sampledata, 
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
     progress = FALSE)
#> Test statistic: 0.06724.             
#> p-value: 0.253.
#> Value of tau used: 0.08311.
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used: 8.
```
Alternatively, the followings are the output when the **full information
approach** is used with the `gurobi` solver to test the hypothesis that
`beta_tgt` is 0.375.
```r
dkqs(df = sampledata, 
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
     progress = FALSE)
#> Test statistic: 0.01746.             
#> p-value: 0.253.
#> Value of tau used: 0.08311.
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used: 8.
```
The results from the two approach should give the same *p*-value while
the test statistic will be different as different moments are considered
in the problem.

The message generated at the end of the function can be re-generated by
`print(dkqs_p)`. All message displayed when the function is calculating
the *p*-value can be re-generated by `summary(dkqs_p)`.

On the other hand, the computational time for using multiple cores
should be shorter than using a single core for a large number of
bootstraps. This is illustrated by the example below:
```r
dkqs_args = list(df = sampledata, 
                 A_obs = A_obs_full, 
                 A_tgt = A_target, 
                 func_obs = func_full_info, 
                 beta_tgt = 0.375, 
                 bs_seed = 1,
                 bs_num = 3000,
                 p_sig = 3,
                 tau_input = tau,
                 solver = "gurobi",
                 cores = 1,
                 progress = FALSE)

# Run dkqs with one core
t10 = Sys.time()
do.call(dkqs, dkqs_args)
#> Test statistic: 0.01746.             
#> p-value: 0.228.
#> Value of tau used: 0.08311.
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used: 1.
t11 = Sys.time()
time1 = t11 - t10

# Run dkqs with eight cores
dkqs_args$cores = 8
t80 = Sys.time()
do.call(dkqs, dkqs_args)
#> Test statistic: 0.01746.             
#> p-value: 0.228.
#> Value of tau used: 0.08311.
#> Linear and quadratic programming solver used: gurobi.
#> Number of cores used: 8.
t81 = Sys.time()
time8 = t81 - t80

# Print the time used
print(sprintf("Time used with 1 core: %s", time1))
#> [1] "Time used with 1 core: 10.4534320831299"
print(sprintf("Time used with 8 cores: %s", time8))
#> [1] "Time used with 8 cores: 4.21064305305481"
```
Test 2: Subsampling
-------------------

The function `subsample` carries out the test using the subsampling
procedure. The syntax is similar to the `dkqs` function. The subsampling
procedure is similar to the bootstrapping procedure except that each
subsample drawn has to have size than the original sample and the draws
are drawn without replacement.

#### Specifying the Functions

Same as the `dkqs` procedure, users need to specify the functions for
the estimator of beta and the variance estimator. The details on how to
construct the function to specify the beta vector can be found <span
id="function">here</span>. In addition, users need to specify the
function to generate the asymptotic variance of the estimator. The order
of the matrix that represent the asymptotic variance has to be the same
as the length of the beta vector.

If the asymptotic variance is an identity matrix, then it can be
specified as follows:
```r
var_full_info <- function(df){
  len = length(unique(df[,"Y"]))
  return(diag(len))
}
```
#### Specifying the Parameters

In this procedure, we need to specify the shape constraints. They can be
defined as follows:
```r
A_shp_full = matrix(rep(1, ncol(A_obs_full)), nrow = 1)
beta_shp = 1
```
#### Syntax

The `subsample` command has the following syntax:
``` r
subsample(df = sampledata, 
          A_obs = A_obs_full, 
          func_obs = func_full_info,
          func_var = var_full_info,
          A_shp = A_shp_full,
          beta_shp = beta_shp,
          A_tgt = A_target,
          beta_tgt = 0.375,
          bs_num = 100,
          bs_seed = 1,
          p_sig = 3,
          solver = "gurobi",
          cores = 8,
          lnorm = 2,
          phi = 2/3,
          progress = FALSE)
```
where

-   `func_var` refers to the function that generates the asymptotic
    variance of observed beta.
-   `A_shp` refers to the matrix for the shape constraints.
-   `beta_shp` refers to the RHS vector for the shape constraints.
-   `lnorm` refers to the norm used in the objective function.
-   `phi` refers to the parameter in determining the size of the
    subsample.

The remaining arguments have the same definition as that has been
defined earlier, which can be found <span id="syntax">here</span>.
Similar to the `dkqs` package, a parallel programming option is
available to increase the computational speed.

#### Sample Output

The following is a sample output of the `subsample` procedure.
```r
subsample(df = sampledata, 
          A_obs = A_obs_full, 
          func_obs = func_full_info,
          func_var = var_full_info,
          A_shp = A_shp_full,
          beta_shp = beta_shp,
          A_tgt = A_target,
          beta_tgt = 0.375,
          bs_num = 100,
          bs_seed = 1,
          p_sig = 3,
          solver = "gurobi",
          cores = 8,
          lnorm = 2,
          phi = 2/3,
          progress = FALSE)
#> Call:
#> subsample(df = sampledata, A_obs = A_obs_full, func_obs = func_full_info, 
#>     func_var = var_full_info, A_shp = A_shp_full, beta_shp = beta_shp, 
#>     A_tgt = A_target, beta_tgt = 0.375, bs_seed = 1, bs_num = 100, 
#>     p_sig = 3, solver = "gurobi", cores = 8, lnorm = 2, phi = 2/3, 
#>     progress = FALSE)
#> 
#> The null hypothesis cannot be rejected at the 5% level.
#> 
#> Test statistic: 0.00055.             
#> p-value: 0.411.
#> Linear and quadratic programming solver used: gurobi.
#> Norm used in the optimization problem: L2-norm.
#> Number of cores used: 8.
```
Constructing Confidence Intervals
---------------------------------

Apart from conducting inference on the tests, the `linearprog` package
can also construct confidence intervals via the `invertci` function. In
this section, I will illustrate the use of the `invertci` function via
the `dkqs` test.

### Constructing a Single Confidence Interval

The command for constructing the two-sided confidence interval for a
statistical test in the `linearprog` packge is called `invertci`. The
confidence interval is constructed by evaluating the *p*-value of a test
and applying the biscetion method.

#### Syntax

The syntax of the `invertci` function is as follows:
``` r
invertci(f = dkqs, 
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

To use the `invertci` function, the arguments for the function of the
test statistic has to be specified and passed to `farg`. For instance,
if the test `dkqs` is used, then the arguments can be defined as
follows:
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
can store them in a data frame and pass it to the function `invertci`.
The requirement for the data frame is as follows:

-   The data frame can only has two columns. The first column is `point`
    (which contains the values of betas that has been evaluated) and the
    second column is `value` (which corresponds to the *p*-values being
    evaluated).
-   The data frame can only contain numeric values.

#### Output

The following shows a sample output of the function `invertci` that is
used to the confidence interval for the test `dkqs` with significance
level 0.05.

``` r
invertci(f = dkqs, 
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
#>  < Constructing confidence interval for alpha = 0.05 >
#> 
#>  === Computing upper bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point  p-value     Reject?
#>  Left end pt.     0.60000        NA            0.60000     0.78000     FALSE  
#>  Right end pt.    NA             1.00000       1.00000     0.00000     TRUE   
#>  1                0.60000        1.00000       0.80000     0.00000     TRUE   
#>  2                0.60000        0.80000       0.70000     0.00000     TRUE   
#>  3                0.60000        0.70000       0.65000     0.05000     FALSE  
#>  4                0.65000        0.70000       0.67500     0.00000     TRUE   
#>  5                0.65000        0.67500       0.66250     0.00000     TRUE   
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  === Computing lower bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point  p-value     Reject?
#>  Left end pt.     0.00000        NA            0.00000     0.00000     TRUE   
#>  Right end pt.    NA             0.40000       0.40000     0.79000     FALSE  
#>  1                0.00000        0.40000       0.20000     0.00000     TRUE   
#>  2                0.20000        0.40000       0.30000     0.00000     TRUE   
#>  3                0.30000        0.40000       0.35000     0.00000     TRUE   
#>  4                0.35000        0.40000       0.37500     0.26000     FALSE  
#>  5                0.35000        0.37500       0.36250     0.03000     FALSE  
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  < Significance level: 0.05 >
#>  Total number of iterations: 10.
#>  Confidence interval: [0.35625, 0.65625].
```

The message generated at the end of the function can be re-generated by
`print(invertci_dkqs)`. All message displayed when the function is
constructing the confidence interval can be re-generated by
`summary(invertci_dkqs)`.

### Constructing Multiple Confidence Intervals

The function `invertci` can also be used to generate multiple confidence
intervals if the argument `alpha` is a vector. For instnace, the
following code produces confidence intervals for alpha equals 0.01, 0.05
and 0.1.

``` r
invertci(f = dkqs, 
         farg = dkqs_farg, 
         alpha = c(0.01, 0.05, 0.1), 
         lb0 = 0, 
         lb1 = 0.4, 
         ub0 = 1, 
         ub1 = 0.6, 
         tol = 0.001, 
         max_iter = 5, 
         df_ci = NULL, 
         progress = TRUE)
#>  < Constructing confidence interval for alpha = 0.01 >
#> 
#>  === Computing upper bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point  p-value     Reject?
#>  Left end pt.     0.60000        NA            0.60000     0.78000     FALSE  
#>  Right end pt.    NA             1.00000       1.00000     0.00000     TRUE   
#>  1                0.60000        1.00000       0.80000     0.00000     TRUE   
#>  2                0.60000        0.80000       0.70000     0.00000     TRUE   
#>  3                0.60000        0.70000       0.65000     0.05000     FALSE  
#>  4                0.65000        0.70000       0.67500     0.00000     TRUE   
#>  5                0.65000        0.67500       0.66250     0.00000     TRUE   
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  === Computing lower bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point  p-value     Reject?
#>  Left end pt.     0.00000        NA            0.00000     0.00000     TRUE   
#>  Right end pt.    NA             0.40000       0.40000     0.79000     FALSE  
#>  1                0.00000        0.40000       0.20000     0.00000     TRUE   
#>  2                0.20000        0.40000       0.30000     0.00000     TRUE   
#>  3                0.30000        0.40000       0.35000     0.00000     TRUE   
#>  4                0.35000        0.40000       0.37500     0.26000     FALSE  
#>  5                0.35000        0.37500       0.36250     0.03000     FALSE  
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  < Constructing confidence interval for alpha = 0.05 >
#> 
#>  === Computing upper bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point  p-value     Reject?
#>  Left end pt.     0.60000        NA            0.60000     0.78000     FALSE  
#>  Right end pt.    NA             1.00000       1.00000     0.00000     TRUE   
#>  1                0.60000        1.00000       0.80000     0.00000     TRUE   
#>  2                0.60000        0.80000       0.70000     0.00000     TRUE   
#>  3                0.60000        0.70000       0.65000     0.05000     FALSE  
#>  4                0.65000        0.70000       0.67500     0.00000     TRUE   
#>  5                0.65000        0.67500       0.66250     0.00000     TRUE   
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  === Computing lower bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point    p-value     Reject?
#>  Left end pt.     0.00000        NA            0.00000       0.00000     TRUE   
#>  Right end pt.    NA             0.40000       0.40000       0.79000     FALSE  
#>  1                0.00000        0.40000       0.20000       0.00000     TRUE   
#>  2                0.20000        0.40000       0.30000       0.00000     TRUE   
#>  3                0.30000        0.40000       0.35000       0.00000     TRUE   
#>  4                0.35000        0.40000       0.37500       0.26000     FALSE  
#>  5                0.35000        0.37500       0.36250       0.03000     FALSE  
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  < Constructing confidence interval for alpha = 0.1 >
#> 
#>  === Computing upper bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point    p-value     Reject?
#>  Left end pt.     0.60000        NA            0.60000       0.78000     FALSE  
#>  Right end pt.    NA             1.00000       1.00000       0.00000     TRUE   
#>  1                0.60000        1.00000       0.80000       0.00000     TRUE   
#>  2                0.60000        0.80000       0.70000       0.00000     TRUE   
#>  3                0.60000        0.70000       0.65000       0.05000     FALSE  
#>  4                0.65000        0.70000       0.67500       0.00000     TRUE   
#>  5                0.65000        0.67500       0.66250       0.00000     TRUE   
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  === Computing lower bound of confidence interval ===
#>  Iteration        Lower bound    Upper bound   Test point    p-value     Reject?
#>  Left end pt.     0.00000        NA            0.00000       0.00000     TRUE   
#>  Right end pt.    NA             0.40000       0.40000       0.79000     FALSE  
#>  1                0.00000        0.40000       0.20000       0.00000     TRUE   
#>  2                0.20000        0.40000       0.30000       0.00000     TRUE   
#>  3                0.30000        0.40000       0.35000       0.00000     TRUE   
#>  4                0.35000        0.40000       0.37500       0.26000     FALSE  
#>  5                0.35000        0.37500       0.36250       0.03000     TRUE   
#>  >>> Reached the maximum number of iterations. Bisection method is completed.
#> 
#>  < Significance level: 0.01 >
#>  Total number of iterations: 10.
#>  Confidence interval: [0.35625, 0.65625].
#>  --------------------------------------
#>  < Significance level: 0.05 >
#>  Total number of iterations: 10.
#>  Confidence interval: [0.35625, 0.65625].
#>  --------------------------------------
#>  < Significance level: 0.1 >
#>  Total number of iterations: 10.
#>  Confidence interval: [0.36875, 0.65625].
```

Constructing Bounds subject to Shape Constraints
------------------------------------------------

Another key application of the `linearprog` package is to construct the
bounds estimates subject to certain shape constraints. This is obtained
by the `estbounds` function. The linear program for obtaining the exact
bounds subject to shape constraints may not be necessarily feasible.
Hence, an `estimate` option for the package is available for the user.
The estimation can be conducted via L1-norm or L2-norm.

#### Syntax

The `estbounds` command has the following syntax:
``` r
estbounds(df = sampledata,
          func_obs = func_full_info,
          A_obs = A_obs_full,
          A_tgt = A_target,
          A_shp_eq = A_shp_eq_dkqs,
          A_shp_ineq = A_shp_ineq_dkqs,
          beta_shp_eq = beta_shp_eq_dkqs,
          beta_shp_ineq = beta_shp_ineq_dkqs,
          kappa = 1e-10,
          lnorm = 2,
          solver = "gurobi",
          estimate = FALSE,
          progress = TRUE)
```
where

-   `A_shp_eq` refers to the matrix representing equality shape
    constraints.
-   `A_shp_ineq` refers to the matrix representing inequality shape
    constraints.
-   `beta_shp_eq` refers to the RHS vector in equality shape
    constraints.
-   `beta_shp_ineq` refers to the RHS vector in inequality shape
    constraints.
-   `lnorm` refers to the norm used in the optimization problem. The
    norms that are supported by this function are L1-norm and L2-norm.
    See <span id="norm">here</span> for more details.
-   `kappa` refers to the parameter used in the second step of the
    two-step procedure for obtaining the solution subje ct to the shape
    constraints.
-   `estimate` refers to the boolean variable that indicate whether the
    estimated problem should be considered.

The remaining arguments have the same definition as that has been
defined earlier, which can be found <span id="syntax">here</span>.

#### Norms

In constructing the estimated bounds, users are free to choose the
L1-norm or the L2-norm. For the estimation with L2-norm, users need to
choose `gurobi` as the solver. For the estimation with L1-norm, users
can choose one of the following packages as the solver:

-   `gurobi`,
-   `limSolve`,
-   `cplexAPI`,
-   `Rcplex`,
-   `lpSolveAPI`.

#### Specifying the Matrices for Shape Constraints

The shape constraints are characterized by equality and inequality
constraints. Each type of the contraints are represented by a matrix and
a vector.

For instance, the equality constraints may be defined as follows:
``` r
A_shp_eq_dkqs = matrix(rep(1, ncol(A_obs_full)), nrow = 1)
beta_shp_eq_dkqs = c(1)
```
If the other set of contraints is empty, set `NULL` to the matrix and
vector.

    A_shp_ineq_dkqs = NULL
    beta_shp_ineq_dkqs = NULL

#### Sample Output

The following is a sample output when the bounds are estimated with a
tolerance level of `1e-20` with L1-norm and the `gurobi` solver.
``` r
estbounds(df = sampledata,
          func_obs = func_full_info,
          A_obs = A_obs_full,
          A_tgt = A_target,
          A_shp_eq = A_shp_eq_dkqs,
          A_shp_ineq = A_shp_ineq_dkqs,
          beta_shp_eq = beta_shp_eq_dkqs,
          beta_shp_ineq = beta_shp_ineq_dkqs,
          kappa = 1e-20,
          lnorm = 2,
          solver = "gurobi",
          estimate = TRUE,
          progress = TRUE)
#> Call:
#> estbounds(df = sampledata, func_obs = func_full_info, A_obs = A_obs_full, 
#>     A_tgt = A_target, A_shp_eq = A_shp_eq_dkqs, A_shp_ineq = A_shp_ineq_dkqs, 
#>     beta_shp_eq = beta_shp_eq_dkqs, beta_shp_ineq = beta_shp_ineq_dkqs, 
#>     kappa = 1e-20, lnorm = 2, solver = "gurobi", estimate = TRUE, 
#>     progress = TRUE)
#> 
#> Norm used in optimization problem: L2-norm 
#> Estimated bounds subject to shape constraints: [0.38316, 0.63344]
```
The following is the sample output when the actual bounds are
constructed without estimation.
```r
estbounds(df = sampledata,
          func_obs = func_full_info,
          A_obs = A_obs_full,
          A_tgt = A_target,
          A_shp_eq = A_shp_eq_dkqs,
          A_shp_ineq = A_shp_ineq_dkqs,
          beta_shp_eq = beta_shp_eq_dkqs,
          beta_shp_ineq = beta_shp_ineq_dkqs,
          kappa = 1e-10,
          lnorm = 2,
          solver = "gurobi",
          estimate = FALSE,
          progress = TRUE)
#> Call:
#> estbounds(df = sampledata, func_obs = func_full_info, A_obs = A_obs_full, 
#>     A_tgt = A_target, A_shp_eq = A_shp_eq_dkqs, A_shp_ineq = A_shp_ineq_dkqs, 
#>     beta_shp_eq = beta_shp_eq_dkqs, beta_shp_ineq = beta_shp_ineq_dkqs, 
#>     kappa = 1e-10, lnorm = 2, solver = "gurobi", estimate = FALSE, 
#>     progress = TRUE)
#> 
#> True bounds subject to shape constraints: [0.3832, 0.6332]
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
