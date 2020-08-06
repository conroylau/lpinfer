---
title: 'lpinfer: An R Package for Inference in Linear Programs'
author: ''
output:
  github_document:
    keep_html: yes
    toc: yes
    df_print: paged
  pdf_document:
    toc: yes
    df_print: paged
  md_document:
    preserve_yaml: yes
    toc: yes
    df_print: paged  
  html_document:
    keep_md: yes
    toc: yes
    df_print: paged
bibliography: refs.bib
vignette: |
  %\VignetteIndexEntry{lpinfer} 
  %\VignetteEncoding{UTF-8} 
  %\VignetteEngine{knitr::rmarkdown}
---

-   [Introduction](#introduction)
-   [Scope of the Vignette](#scope-of-the-vignette)
-   [Installation and Requirements](#installation-and-requirements)
-   [Sample data](#sample-data)
-   [General Syntax](#general-syntax)
-   [Tests](#tests)
    -   [Test 1: DKQS cone-tightening
        procedure](#test-1-dkqs-cone-tightening-procedure)
    -   [Test 2: Subsampling procedure](#test-2-subsampling-procedure)
    -   [Test 3: FSST procedure](#test-3-fsst-procedure)
-   [Constructing Bounds subject to Shape
    Constraints](#constructing-bounds-subject-to-shape-constraints)
-   [Constructing Confidence
    Intervals](#constructing-confidence-intervals)
-   [Parallel Programming](#parallel-programming)
-   [Help, Feature Requests and Bug
    Reports](#help-feature-requests-and-bug-reports)
-   [References](#references)

Introduction
------------

This package provides a set of methods to conduct inference on
econometrics problems that can be studied by linear programs. Currently,
this package supports the following tests:

1.  Cone-tightening procedure by Deb et al. (2018)
2.  Subsampling procedure
3.  FSST procedure by Fang et al. (2020)

Apart from computing the *p*-values based on the above tests, this
package can also construct confidence intervals and estimate the bounds
of the estimators in linear programs subject to certain shape
constraints.

Scope of the Vignette
---------------------

This vignette is intended as a guide to use the `lpinfer` package
without further explanation for the methods. Readers may refer to each
of the sections for the tests for references to more details about the
methods.

Installation and Requirements
-----------------------------

`lpinfer` can be installed from our GitHub repository via

```r
devtools::install_github("conroylau/lpinfer")
```

To use most of the functions in this `lpinfer` package, one of the
following packages for solving the linear and quadratic programs is
required. There are four options for the solver:

1.  Gurobi and the R package `gurobi` — Gurobi can be downloaded from
    [Gurobi Optimization](https://www.gurobi.com/). A Gurobi software
    license is required. The license can be obtained at no cost for
    academic researchers. The instructions for installing `gurobi` on R
    can be found
    [here](https://cran.r-project.org/web/packages/prioritizr/vignettes/gurobi_installation.html#r-package-installation).

2.  IBM ILOG CPLEX Optimization Studio (CPLEX) and one of the R packages
    below — CPLEX can be downloaded from
    [IBM](https://www.ibm.com/analytics/cplex-optimizer). A CPLEX
    software license is required, which can be obtained at no cost for
    academic researchers. There are two open-source and free R packages
    that uses CPLEX, and users are free to choose one of them. In
    addition, both packages have to be installed on the command line to
    link the package to the correct CPLEX library. The two packages’
    name and installation instructions are as follows:

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

Sample data
-----------

The classical missing data problem due to Manski (1989) is used as a
running example throughout this vignette to demonstrate the commands in
this `lpinfer` package. We consider this problem here because the sharp
bounds for the identified set of the expected values can be constructed
by linear programs.

In this package, a sample simulated data set on the missing data problem
is included. This can be obtained by `sampledata`. This data set
contains 1,000 observations with 2 columns. The following shows the
first 10 observations of the simulated data set:

```r
library(lpinfer)
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
    from 0 to 1 with increment 0.1.
-   `D` is a binary treatment variable where *Y*<sub>*i*</sub> is
    observed for *D*<sub>*i*</sub> = 1 and not observed for
    *D*<sub>*i*</sub> = 0.

General Syntax
--------------

In general, each of the tests in the `lpinfer` requires a data set, a
`lpmodel` object and some tuning parameters. The tuning parameters are
different for each of the test and will be explained in later sections.
For the `lpmodel` object, it consists of five components:

-   `A.obs`
-   `A.shp`
-   `A.tgt`
-   `beta.obs`
-   `beta.shp`

Here, `A.obs` and `beta.obs` refers to the matrix of coefficients and
the RHS variables for the “observed constraints” in the linear program.
`A.shp` and `beta.shp` refers to the matrix of coefficients and the RHS
variables for the shape constraints. `A.tgt` refers to the matrix of
coefficients for the linear constraints that we wish to test with
`beta.tgt`. Since the object `beta.tgt` is a parameter that is being
tested, it will be defined outside the `lpmodel` object.

Note that not all procedures will require all five components to be
present in the `lpmodel` object. For instance, `A.shp` and `beta.shp`
are not required in the `dkqs` procedure.

#### Standard form

For consistency, the tests in this `lpinfer` package assumes that all of
the components in the `lpmodel` object represent constraints that are  
in standard form. Hence, the matrices in the `lpmodel` object are
representing equality constraints. To impose inequality constraints,
users need to first convert them to standard form by adding appropriate
slack and surplus variables.

This can be done easily with the `standard.lpmodel` function in this
`lpinfer` package with an object in the `lpmodel.natural` class. The
`lpmodel.natural` class consists of the following eight components:

-   `A.obs`
-   `A.shp`
-   `A.tgt`
-   `beta.obs`
-   `beta.shp`
-   `sense.shp`
-   `x.lb`
-   `x.ub`

The first five components are the same as the `lpmodel` class. The
`sense.shp` vector stores the sense of the shape constraints (i.e. `>=`,
`=` or `<=`). The two components `x.lb` and `x.ub` correspond to the
lower and upper bound of the `x` vector respectively.

The procedure to convert an object from the `lpmodel.natural` class into
the `lpmodel` class is as follows:

1.  Create an object in the `lpmodel.natural` class.
2.  Apply the `standard.lpmodel` function on this object.

Below is an example of using the `standard.lpmodel` function on the
object in the `lpmodel.natural` class, which is denoted by `lpmn0` in
this example.

```r
### Step 1: Create an object in the `lpmodel.natural` class
# Obs
Aobs0 <- matrix(c(1, 2), nrow = 1)
bobs0 <- c(10)

# Shp
Ashp0 <- matrix(c(3, 4, 5, 6), nrow = 2, byrow = TRUE)
bshp0 <- matrix(c(15, 100))
sshp0 <- matrix(c(">=", "<="))

# Tgt
Atgt0 <- matrix(c(1, 1), nrow = 1)

# Upper bounds
xub0 <- c(200, 200)

# Lower bounds
xlb0 <- c(0.1, 0.1)

# Formulate the `lpmodel.natural` object
lpmn0 <- lpmodel.natural(A.obs = Aobs0,
                         A.shp = Ashp0,
                         A.tgt = Atgt0,
                         beta.obs = bobs0,
                         beta.shp = bshp0,
                         sense.shp = sshp0,
                         x.ub = xub0,
                         x.lb = xlb0)
```

Then, these constraints can be transformed into a `lpmodel` object in
standard form as follows:

```r
### Step 2: Apply the `standard.lpmodel` function
lpm1 <- standard.lpmodel(lpmn0)
```

The updated `A.shp` and `beta.shp` matrices reflect the updated shape
constraints and the lower and upper bounds in standard form by
incorporating slack and surplus variables. The two matrices can be
viewed as follows:

```r
print(lpm1$A.shp)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
#> [1,]    5    6    1    0    0    0    0    0
#> [2,]    1    0    0    1    0    0    0    0
#> [3,]    0    1    0    0    1    0    0    0
#> [4,]    3    4    0    0    0   -1    0    0
#> [5,]    1    0    0    0    0    0   -1    0
#> [6,]    0    1    0    0    0    0    0   -1
print(lpm1$beta.shp)
#>       [,1]
#> [1,] 100.0
#> [2,] 200.0
#> [3,] 200.0
#> [4,]  15.0
#> [5,]   0.1
#> [6,]   0.1
```

-   The first row corresponds to the `<=` constraint in the original
    `Ashp0` and `bshp0` matrix.
-   The second and third rows correspond to the upper bounds.
-   The fourth row corresponds to the `>=` constraint in the original
    `Ashp0` and `bshp0` matrix.
-   The last two rows correspond to the lower bounds.

The resulting object `lpm1` is already in standard form and can be
passed to the tests in this `lpinfer` package.

#### Deterministic or stochastic components

Depending on the testing procedure, the five components in the `lpmodel`
object can be either deterministic or stochastic. The component is said
to be deterministic if it remains unchanged in the bootstrap procedure.
The component is said to be stochastic if it changes (and depends on the
bootstrap data) in the bootstrap procedure.

-   If the component is deterministic, then the component is typically a
    `matrix` or a `data.frame`.

-   If the component is stochastic, then it is represented by a
    `function` or a `list`. If the component is a `function`, then it
    will be re-evaluated based on the bootstrap data in each bootstrap
    draw. If the component is a `list`, then each object in the list
    represents the component in the bootstrap draw.

#### Remarks for the component being a function

If the component is a `function`, then it has to fulfill the following
requirements:

-   The function’s only argument is the data set. The function needs to
    accept data sets in the class `data.frame`.
-   The function can return either one or two objects. If it returns one
    object, then it has to be a vector or a column matrix that
    represents the estimator for `beta.obs`. If it has two objects, then
    it returns the estimator for `beta.obs` and a square matrix that
    refers to the estimator of the asymptotic variance of `beta.obs`.

#### Demonstration

The following is an example on how to construct each of the components
in the `lpmodel` object and how to assign the `lpmodel` object based on
the example data mentioned [here](#sample-data):

```r
# Extract relevant information from data
N <- nrow(sampledata)
J <- length(unique(sampledata[,"Y"])) - 1
J1 <- J + 1

# Construct A.obs
Aobs.full <- cbind(matrix(rep(0, J1 * J1), nrow = J1), diag(1, J1))

# Construct A.tgt
yp <- seq(0, 1, 1/J)
Atgt <- matrix(c(yp, yp), nrow = 1)

# Construct A.shp
Ashp <- matrix(rep(1, ncol(Aobs.full)), nrow = 1)

# Construct beta.obs (a function)
betaobs.fullinfo <- function(data){
  beta <- NULL
  y.list <- sort(unique(data[,"Y"]))
  n <- dim(data)[1]
  yn <- length(y.list)
  for (i in 1:yn) {
    beta.i <- sum((data[,"Y"] == y.list[i]) * (data[,"D"] == 1))/n
    beta <- c(beta, c(beta.i))
  }
  beta <- as.matrix(beta)
  return(list(beta = beta,
              var = diag(yn)))
}

# Define the lpmodel object
lpm.full <- lpmodel(A.obs    = Aobs.full,
                    A.tgt    = Atgt,
                    A.shp    = Ashp,
                    beta.obs = betaobs.fullinfo,
                    beta.shp = c(1))
```

In the above, the `betaobs.fullinfo` function returns two information:
`beta` is a vector of length that is equal to the number of distinct
observations for *Y* where element *i* of the vector refers to the
probability that the corresponding value of *y*<sub>*i*</sub> is
observed. The `var` object is the estimator of the asymptotic variance
of `beta`, which is assumed to be an identity matrix here for
illustration purpose.

Tests
-----

### Test 1: DKQS cone-tightening procedure

The `dkqs` function in this `lpinfer` package is used to conduct the
cone-tightening procedure that is proposed by Deb et al. (2018). For
details of the procedure, readers may refer to section 4.2 and the
appendix of Deb et al. (2018).

#### Syntax

The `dkqs` command has the following syntax:

```r
dkqs(data = sampledata,
     lpmodel = lpm.full,
     beta.tgt = .375,
     R = 100,
     tau = sqrt(log(N)/N),
     solver = "gurobi",
     progress = FALSE)
```

where

-   `data` refers to the data set.
-   `lpmodel` refers to the `lpmodel` object.
-   `beta.tgt` refers to the parameter that is being tested.
-   `R` refers to the total number of bootstraps.
-   `tau` refers to the tuning parameter tau. This will be explained
    [here](#choosing-the-tau-parameter). It can be a scalar or a vector.
-   `progress` refers to the boolean variable for whether the progress
    bar should be printed in the testing procedure.

#### Choosing the tau parameter

Depending on the value of `tau`, the main quadratic program in Deb et
al. (2018) is not always feasible. The details on how to choose the
maximum tau such that the `dkqs` procedure is feasible can be found in
the supplemental appendix of Kamat (2019).

In this module, we follow the procedure by Kamat (2019) to pick the
largest feasible `tau`. If the value of `tau` chosen by the user is
infeasible, then the module will evaluate the procedure using the
largest possible `tau`. Otherwise, the module will evaluate the
procedure with the given `tau`.

If `tau` is not specified by the user, the procedure will directly use
the maximum feasible `tau`.

On the other hand, this package provides the flexibility to the users to
pass multiple tuning parameters at the same time. Hence, `tau` can be a
vector as well.

#### Components in `lpmodel`

The following table summarizes whether the components in `lpmodel` can
be deterministic or stochastic:

<table>
<thead>
<tr class="header">
<th>Component in <code>lpmodel</code></th>
<th style="text-align: left;">Property</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>A.obs</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
<tr class="even">
<td><code>A.shp</code></td>
<td style="text-align: left;">Not used</td>
</tr>
<tr class="odd">
<td><code>A.tgt</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
<tr class="even">
<td><code>beta.obs</code></td>
<td style="text-align: left;">Stochastic</td>
</tr>
<tr class="odd">
<td><code>beta.shp</code></td>
<td style="text-align: left;">Not used</td>
</tr>
</tbody>
</table>

`A.shp` and `beta.shp` will be ignored in the `dkqs` procedure.

#### Data structure

The function defined in this function needs to generate an estimator for
the `beta.obs` object.

In addition, this procedure has an additional requirement on the data
structure: the column that represents the outcome variable has to be
called either “Y” or “y”.

#### Example

The following is what happens when the above code is run:

```r
set.seed(1)
dkqs.full1 <- dkqs(data = sampledata,
                   lpmodel = lpm.full,
                   beta.tgt = .375,
                   R = 100,
                   tau = sqrt(log(N)/N),
                   solver = "gurobi",
                   progress = FALSE)
print(dkqs.full1)
#>  p-value: 0.22
```

As noted in the syntax section, we provide the flexibility for users to
conduct inference with multiple tuning parameters at the same time. The
following is an example when the user specifies multiple taus in the
`dkqs` procedure:

```r
set.seed(1)
dkqs.full2 <- dkqs(data = sampledata,
                   lpmodel = lpm.full,
                   beta.tgt = .375,
                   R = 100,
                   tau = c(.1, .2, .3, .5),
                   solver = "gurobi",
                   progress = FALSE)
print(dkqs.full2)
#>  tau p-value
#>  0.1    0.22
#>  0.2    0.22
#>  0.3    0.22
#>  0.5    0.29
```

Users can get a more detailed summary of the results when applying the
`summary` command on the resulting object:

```r
summary(dkqs.full2)
#>  tau p-value
#>  0.1    0.22
#>  0.2    0.22
#>  0.3    0.22
#>  0.5    0.29
#>  Maximum feasible tau: 0.76923
#>  Test statistic: 0.01746
#>  Solver used: gurobi
#>  Number of successful bootstrap replications: 100
```

#### Alternative approach

The approach that has been demonstrated so far is known as the full
information approach. The inference procedures in this `lpinfer` package
is not limited to any single approach. In this section, we demonstrate
how to use the two moments approach. This alternative method is known as
the two moments approach because the two moments
**E**\[*Y*<sub>*i*</sub>\] and
**E**\[*Y*<sub>*i*</sub>*D*<sub>*i*</sub>\] are used in the inference
procedure.

Similar to what has been demonstrated in the [example](#demonstration)
earlier, some of the components has to be updated in the `lpmodel`
object. This can be done as follows:

```r
# Construct A.obs
Aobs.twom <- matrix(c(rep(0,J1), yp, rep(0,J1), rep(1, J1)), nrow = 2,
                     byrow = TRUE)

# Construct beta.obs (a function)
betaobs.twom <- function(data){
  beta <- matrix(c(0,0), nrow = 2)
  n <- dim(data)[1]
  beta[1] <- sum(data[,"Y"] * data[,"D"])/n
  beta[2] <- sum(data[,"D"])/n
  return(list(beta = beta,
              var = diag(2)))
}

# Define the lpmodel object
lpm.twom <- lpmodel(A.obs    = Aobs.twom,
                    A.tgt    = Atgt,
                    A.shp    = Ashp,
                    beta.obs = betaobs.twom,
                    beta.shp = c(1))
```

The `dkqs` procedure can be ran in the same fashion as before. The
following is an example on how this can be ran:

```r
set.seed(1)
dkqs.full3 <- dkqs(data = sampledata,
                   lpmodel = lpm.twom,
                   beta.tgt = .375,
                   R = 100,
                   tau =  sqrt(log(N)/N),
                   solver = "gurobi",
                   progress = FALSE)
print(dkqs.full3)
#>  p-value: 0.22
```

As shown above, we get the same *p*-value as before.

### Test 2: Subsampling procedure

The `subsample` function in this package carries out the test using the
subsampling procedure.

#### Syntax

The `subsample` command has the following syntax:

```r
subsample(data = sampledata, 
          lpmodel = lpm.full,
          beta.tgt = 0.375,
          R = 100,
          solver = "gurobi",
          norm = 2,
          phi = 2/3,
          replace = FALSE,
          progress = FALSE)
```

where

-   `phi` refers to the parameter that controls the size of each
    subsample. This will be further explained
    [here](#choosing-the-phi-and-replace-parameter).
-   `replace` refers to the boolean variable to indicate whether the
    function samples the data with or without replacement. This will be
    further explained [here](#choosing-the-phi-and-replace-parameter).
-   `norm` refers to the norm used in the objective function.

The rest of the arguments are the same as that in the `dkqs` procedure.

#### Choosing the `phi` and `replace` parameter

The `phi` parameter is a parameter to control the size of each
subsample. The sample size of the original data to the power `phi` is
the size of each subsample. On the other hand, the `replace` parameter
is used to indicate whether the function samples the data with or
without replacement.

<table>
<thead>
<tr class="header">
<th><code>replace</code></th>
<th style="text-align: left;"><code>phi</code></th>
<th style="text-align: left;">Meaning</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>FALSE</code></td>
<td style="text-align: left;">Any number in the interval (0, 1)</td>
<td style="text-align: left;">Subsample</td>
</tr>
<tr class="even">
<td><code>TRUE</code></td>
<td style="text-align: left;">Equals to 1</td>
<td style="text-align: left;">Bootstrap</td>
</tr>
<tr class="odd">
<td><code>TRUE</code></td>
<td style="text-align: left;">Any number in the interval (0, 1)</td>
<td style="text-align: left;"><em>m</em> out of <em>n</em> bootstrap</td>
</tr>
</tbody>
</table>

Note that users cannot specify `phi` as 1 when `replace` is set to
`FALSE` because it will be generating the exactly same set of data in
every subsample draw.

#### Components in `lpmodel`

The following table summarizes whether the components in `lpmodel` can
be deterministic or stochastic:

<table>
<thead>
<tr class="header">
<th>Component in <code>lpmodel</code></th>
<th style="text-align: left;">Property</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>A.obs</code></td>
<td style="text-align: left;">Deterministic or stochastic</td>
</tr>
<tr class="even">
<td><code>A.shp</code></td>
<td style="text-align: left;">Deterministic or stochastic</td>
</tr>
<tr class="odd">
<td><code>A.tgt</code></td>
<td style="text-align: left;">Deterministic or stochastic</td>
</tr>
<tr class="even">
<td><code>beta.obs</code></td>
<td style="text-align: left;">Stochastic</td>
</tr>
<tr class="odd">
<td><code>beta.shp</code></td>
<td style="text-align: left;">Deterministic or stochastic</td>
</tr>
</tbody>
</table>

#### Example

The following is what happens when the above code is run:

```r
set.seed(1)
subsample.full <- subsample(data = sampledata, 
                            lpmodel = lpm.full,
                            beta.tgt = 0.375,
                            R = 100,
                            solver = "gurobi",
                            norm = 2,
                            phi = 2/3,
                            replace = FALSE,
                            progress = FALSE)
print(subsample.full)
#> p-value: 0.3
```

Again, more detailed information can be extracted via the `summary`
command:

```r
summary(subsample.full)
#> p-value: 0.3
#> Test statistic: 0.01746
#> Solver used: gurobi
#> Norm used: 2
#> Phi used: 0.66667
#> Size of each subsample: 99
#> Number of successful bootstrap replications: 100
```

As indicated [earlier](#choosing-the-phi-and-replace-parameter), the
`subsample` command can perform the bootstrap and *m* out of *n*
boostrap procedures. They are illustrated as follows.

The following is an exmaple of performing a bootstrapping procedure:

```r
set.seed(1)
subsample.bootstrap <- subsample(data = sampledata, 
                                 lpmodel = lpm.full,
                                 beta.tgt = 0.375,
                                 R = 100,
                                 solver = "gurobi",
                                 norm = 2,
                                 phi = 1,
                                 replace = TRUE,
                                 progress = FALSE)
print(subsample.bootstrap)
#> p-value: 0.48
```

The following is an exmaple of performing a *m* out of *n* bootstrapping
procedure:

```r
set.seed(1)
subsample.bootstrap2 <- subsample(data = sampledata, 
                                  lpmodel = lpm.full,
                                  beta.tgt = 0.375,
                                  R = 100,
                                  solver = "gurobi",
                                  norm = 2,
                                  phi = 2/3,
                                  replace = TRUE,
                                  progress = FALSE)
print(subsample.bootstrap2)
#> p-value: 0.3
```

### Test 3: FSST procedure

The `fsst` function in this `lpinfer` package is used to conduct the
testing procedure by Fang et al. (2020).

#### Syntax

The `fsst` command has the following syntax:

```r
fsst(data = sampledata, 
     lpmodel = lpm.full,
     beta.tgt = 0.375,
     R = 100,
     lambda = 0.5,
     rho = 1e-4,
     n = nrow(sampledata),
     weight.matrix = "diag",
     solver = "gurobi",
     progress = FALSE)
```

where

-   `lambda` refers to the tuning parameter that affects test statistics
    in the bootstrap cone components. Users can pass multiple `lambda`s
    in this argument.
-   `rho` refers to the parameter used to studentize the variance
    matrices in the FSST procedure.
-   `n` is optional if `data` is passed. `n` is a variable that refers
    to the number of rows of `data`. If the `beta.obs` in `lpmodel` is a
    `list`, then users can skip `data` and pass the number of rows of
    the `data` as `n` instead.
-   `weight.matrix` is a string that determines the weighting matrix.
    The details can be found [here](#weight-matrix).

The rest of the arguments are the same as that in the `dkqs` procedure.

#### Weighting matrix

This procedure provides three options for the weighting matrix in the
FSST procedure:

<table>
<thead>
<tr class="header">
<th><code>weight.matrix</code></th>
<th style="text-align: left;">Definition</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>identity</code></td>
<td style="text-align: left;">Identity matrix</td>
</tr>
<tr class="even">
<td><code>avar</code></td>
<td style="text-align: left;">Inverse of the asymptotic variance of <code>beta.obs</code></td>
</tr>
<tr class="odd">
<td><code>diag</code></td>
<td style="text-align: left;">Diagonal matrix with elements equal to the diagonal entries of <code>avar</code></td>
</tr>
</tbody>
</table>

#### Components in `lpmodel`

The following table summarizes whether the components in `lpmodel` can
be deterministic or stochastic:

<table>
<thead>
<tr class="header">
<th>Component in <code>lpmodel</code></th>
<th style="text-align: left;">Property</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>A.obs</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
<tr class="even">
<td><code>A.shp</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
<tr class="odd">
<td><code>A.tgt</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
<tr class="even">
<td><code>beta.obs</code></td>
<td style="text-align: left;">Stochastic</td>
</tr>
<tr class="odd">
<td><code>beta.shp</code></td>
<td style="text-align: left;">Deterministic</td>
</tr>
</tbody>
</table>

#### Example

The following is what happens when the above code is run:

```r
set.seed(1)
fsst.full1 <- fsst(data = sampledata, 
                   lpmodel = lpm.full,
                   beta.tgt = 0.375,
                   R = 100,
                   lambda = 0.5,
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "diag",
                   solver = "gurobi",
                   progress = FALSE)
print(fsst.full1)
#> p-value: 0.22
```

As noted above, the `fsst` procedure provides the flexibility for users
to pass multiple tuning parameters.

```r
set.seed(1)
fsst.full2 <- fsst(data = sampledata, 
                   lpmodel = lpm.full,
                   beta.tgt = 0.375,
                   R = 100,
                   lambda = c(0.1, 0.2, 0.5),
                   rho = 1e-4,
                   n = nrow(sampledata),
                   weight.matrix = "diag",
                   solver = "gurobi",
                   progress = FALSE)
print(fsst.full2)
#>      lambda  p-value
#>      0.1 0.95
#>      0.2 0.51
#>      0.5 0.22
```

Again, more detailed information can be extracted via the `summary`
command:

```r
summary(fsst.full2)
#> 
#> Sample and quantiles of bootstrap test statistics: 
#>                               lambda       0.1     0.2     0.5
#>     Test statistic            Sample   0.04698 0.04698 0.04698
#>                     Bootstrap 99% CV   0.55484 0.31925 0.11116
#>                     Bootstrap 95% CV   0.43154 0.23072  0.0974
#>                     Bootstrap 90% CV   0.35565 0.16116  0.0762
#>               Cone            Sample   0.04698 0.04698 0.04698
#>                     Bootstrap 99% CV   0.55484 0.31925 0.11116
#>                     Bootstrap 95% CV   0.43154 0.23072  0.0974
#>                     Bootstrap 90% CV   0.35565 0.16116  0.0762
#>              Range            Sample   0.00000                
#>                     Bootstrap 99% CV   0.00000                
#>                     Bootstrap 95% CV   0.00000                
#>                     Bootstrap 90% CV   0.00000                
#> 
#> p-values:
#>      lambda      0.1  0.2  0.5
#>      p-value    0.95 0.51 0.22
#> 
#> Solver used: gurobi
#> 
#> Regularization parameters: 
#>    - Input value of rho: 1e-04
#>    - Regularization parameter for the Cone studentization matrix: 0.00033
#> 
#> The asymptotic variance of the observed component of the beta vector is approximated from the function.
#> 
#> Number of successful bootstrap replications: 100
```

Constructing Bounds subject to Shape Constraints
------------------------------------------------

Another key application of the `lpinfer` package is to construct the
bounds estimates subject to certain shape constraints. This is obtained
by the `estbounds` function. The linear program for obtaining the exact
bounds subject to shape constraints may not be necessarily feasible.
Hence, an `estimate` option is available. The estimation can be
conducted via L1-norm or L2-norm.

#### Syntax

```r
estbounds(data = sampledata,
          lpmodel = lpm.full,
          kappa = 1e-5,
          norm = 2,
          solver = "gurobi",
          estimate = TRUE,
          progress = FALSE)
```

where

-   `norm` refers to the norm used in the optimization problem. The
    norms that are supported by this function are L1-norm and L2-norm.
    See [here](#norms) for more details.
-   `kappa` refers to the parameter used in the second step of the
    two-step procedure for obtaining the solution subject to the shape
    constraints.
-   `estimate` refers to the boolean variable that indicate whether the
    estimated problem should be considered.

The rest of the arguments are the same as that in the `dkqs` procedure.

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

#### Example

The following is what happens when the above code is run:

```r
set.seed(1)
estbounds.full <- estbounds(data = sampledata,
                            lpmodel = lpm.full,
                            kappa = 1e-5,
                            norm = 2,
                            solver = "gurobi",
                            estimate = TRUE,
                            progress = FALSE)
print(estbounds.full)
#> Estimated bounds: [0.38316, 0.63344]
```

Again, more detailed information can be extracted via the `summary`
command:

```r
summary(estbounds.full)
#> Estimated bounds: [0.38316, 0.63344] 
#> Norm used: 2 
#> Solver: gurobi
```

Constructing Confidence Intervals
---------------------------------

Apart from conducting inference and estimating the bounds, the `lpinfer`
package can also construct confidence intervals via the `invertci`
function. The confidence interval is constructed by evaluating the
*p*-value of a test and applying the bisection method.

#### Syntax

The syntax of the `invertci` function is as follows:

```r
invertci(f = subsample, 
         farg = subsample.args, 
         alpha = 0.05, 
         lb0 = NULL, 
         lb1 = NULL, 
         ub0 = NULL, 
         ub1 = NULL, 
         tol = 0.0001, 
         max.iter = 20, 
         df_ci = NULL,
         progress = FALSE)
```

where

-   `f` refers to the function that represents a testing procedure.
-   `farg` refers to the list of arguments to be passed to the function
    of testing procedure. The details can be found
    [here](#specifying-the-argument).
-   `alpha` refers to the significance level of the test. Please refer
    to the details [here](#constructing-multiple-confidence-intervals).
-   `lb0` refers to the logical lower bound for the confidence interval.
-   `lb1` refers to the maximum possible lower bound for the confidence
    interval.
-   `ub0` refers to the logical upper bound for the confidence interval.
-   `ub1` refers to the minimum possible upper bound for the confidence
    interval.
-   `tol` refers to the tolerance level in the bisection method.
-   `max.iter` refers to the maximum number of iterations in the
    bisection method.
-   `df_ci` refers to `data.frame` that consists of the points and the
    corresponding *p*-values that have been tested in constructing the
    confidence intervals. The details can be found
    [here](#specifying-the-data-frame-df_ci).
-   `progress` refers to the boolean variable for whether the result
    messages should be displayed in the procedure of constructing
    confidence interval.

#### Specifying the argument

To use the `invertci` function, the arguments for the test statistic has
to be specified and passed to the `farg` argument. For instance, if the
testing procedure `subsample` is used, the arguments can be defined as
follows:

```r
subsample.args <- list(data = sampledata,
                       lpmodel = lpm.full,
                       R = 100,
                       phi = .75,
                       solver = "gurobi",
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

#### Example

The following shows a sample output of the function `invertci` that is
used to the confidence interval for the test `dkqs` with significance
level 0.05.

```r
set.seed(1)
invertci.subsample1 <- invertci(f = subsample,
                                farg = subsample.args, 
                                alpha = 0.05, 
                                lb0 = 0, 
                                lb1 = 0.4, 
                                ub0 = 1, 
                                ub1 = 0.6, 
                                tol = 0.001, 
                                max.iter = 5, 
                                df_ci = NULL, 
                                progress = FALSE)
print(invertci.subsample1)
#> 
#> Confidence interval: [0.35625, 0.66875]
```

The details for each iteration can be obtained in real-time by setting
`progress` as `TRUE` or applying the `summary` command on the resulting
object:

```r
summary(invertci.subsample1)
#> 
#> Significance level: 0.05
#> Confidence interval: [0.35625, 0.66875]
#> 
#> Maximum number of iterations: 5
#> Tolerance level: 0.001
#> 
#> Details:
#> 
#> === Iterations in constructing upper bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.60000            NA      0.60000   0.57000     FALSE
#> Right end pt.            NA       1.00000      1.00000   0.00000      TRUE
#> 1                   0.60000       1.00000      0.80000   0.00000      TRUE
#> 2                   0.60000       0.80000      0.70000   0.00000      TRUE
#> 3                   0.60000       0.70000      0.65000   0.12000     FALSE
#> 4                   0.65000       0.70000      0.67500   0.00000      TRUE
#> 5                   0.65000       0.67500      0.66250   0.03000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> === Iterations in constructing lower bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.00000            NA      0.00000   0.00000      TRUE
#> Right end pt.            NA       0.40000      0.40000   0.46000     FALSE
#> 1                   0.00000       0.40000      0.20000   0.00000      TRUE
#> 2                   0.20000       0.40000      0.30000   0.00000      TRUE
#> 3                   0.30000       0.40000      0.35000   0.02000      TRUE
#> 4                   0.35000       0.40000      0.37500   0.28000     FALSE
#> 5                   0.35000       0.37500      0.36250   0.10000     FALSE
#> Reason for termination: Reached maximum number of iterations
```

#### Constructing Multiple Confidence Intervals

The function `invertci` can also be used to generate multiple confidence
intervals if the argument `alpha` is a vector. For instance, the
following code produces confidence intervals for alpha equals 0.01, 0.05
and 0.1.

```r
set.seed(1)
invertci.subsample2 <- invertci(f = subsample, 
                                farg = subsample.args, 
                                alpha = c(0.01, 0.05, 0.1), 
                                lb0 = 0, 
                                lb1 = 0.4, 
                                ub0 = 1, 
                                ub1 = 0.6, 
                                tol = 0.001, 
                                max.iter = 5, 
                                df_ci = NULL, 
                                progress = FALSE)
print(invertci.subsample2)
#>                                       
#> Significance level Confidence interval
#> 0.01                [0.33125, 0.66875]
#> 0.05                [0.34375, 0.66875]
#> 0.1                 [0.35625, 0.66875]
```

Again, the detailed steps in constructing the confidence intervals can
be obtained as follows:

```r
summary(invertci.subsample2)
#> Maximum number of iterations: 5
#> Tolerance level: 0.001
#>                                       
#> Significance level Confidence interval
#> 0.01                [0.33125, 0.66875]
#> 0.05                [0.34375, 0.66875]
#> 0.1                 [0.35625, 0.66875]
#> 
#> Details:
#> 
#> <Confidence interval for significance level = 0.01>
#> 
#> === Iterations in constructing upper bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.60000            NA      0.60000   0.57000     FALSE
#> Right end pt.            NA       1.00000      1.00000   0.00000      TRUE
#> 1                   0.60000       1.00000      0.80000   0.00000      TRUE
#> 2                   0.60000       0.80000      0.70000   0.00000      TRUE
#> 3                   0.60000       0.70000      0.65000   0.12000     FALSE
#> 4                   0.65000       0.70000      0.67500   0.00000      TRUE
#> 5                   0.65000       0.67500      0.66250   0.03000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> === Iterations in constructing lower bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.00000            NA      0.00000   0.00000      TRUE
#> Right end pt.            NA       0.40000      0.40000   0.46000     FALSE
#> 1                   0.00000       0.40000      0.20000   0.00000      TRUE
#> 2                   0.20000       0.40000      0.30000   0.00000      TRUE
#> 3                   0.30000       0.40000      0.35000   0.02000     FALSE
#> 4                   0.30000       0.35000      0.32500   0.00000      TRUE
#> 5                   0.32500       0.35000      0.33750   0.01000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> <Confidence interval for significance level = 0.05>
#> 
#> === Iterations in constructing upper bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.60000            NA      0.60000   0.48000     FALSE
#> Right end pt.            NA       1.00000      1.00000   0.00000      TRUE
#> 1                   0.60000       1.00000      0.80000   0.00000      TRUE
#> 2                   0.60000       0.80000      0.70000   0.00000      TRUE
#> 3                   0.60000       0.70000      0.65000   0.15000     FALSE
#> 4                   0.65000       0.70000      0.67500   0.01000      TRUE
#> 5                   0.65000       0.67500      0.66250   0.05000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> === Iterations in constructing lower bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.00000            NA      0.00000   0.00000      TRUE
#> Right end pt.            NA       0.40000      0.40000   0.46000     FALSE
#> 1                   0.00000       0.40000      0.20000   0.00000      TRUE
#> 2                   0.20000       0.40000      0.30000   0.00000      TRUE
#> 3                   0.30000       0.40000      0.35000   0.05000     FALSE
#> 4                   0.30000       0.35000      0.32500   0.00000      TRUE
#> 5                   0.32500       0.35000      0.33750   0.00000      TRUE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> <Confidence interval for significance level = 0.1>
#> 
#> === Iterations in constructing upper bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.60000            NA      0.60000   0.55000     FALSE
#> Right end pt.            NA       1.00000      1.00000   0.00000      TRUE
#> 1                   0.60000       1.00000      0.80000   0.00000      TRUE
#> 2                   0.60000       0.80000      0.70000   0.00000      TRUE
#> 3                   0.60000       0.70000      0.65000   0.10000     FALSE
#> 4                   0.65000       0.70000      0.67500   0.02000      TRUE
#> 5                   0.65000       0.67500      0.66250   0.05000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> === Iterations in constructing lower bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.00000            NA      0.00000   0.00000      TRUE
#> Right end pt.            NA       0.40000      0.40000   0.44000     FALSE
#> 1                   0.00000       0.40000      0.20000   0.00000      TRUE
#> 2                   0.20000       0.40000      0.30000   0.00000      TRUE
#> 3                   0.30000       0.40000      0.35000   0.02000      TRUE
#> 4                   0.35000       0.40000      0.37500   0.27000     FALSE
#> 5                   0.35000       0.37500      0.36250   0.10000     FALSE
#> Reason for termination: Reached maximum number of iterations
```

Since the details of the iterations could potentially be a long list of
outputs, users may specify a particular list of output that they wish to
read in the `summary` command. For instance, to only print the details
of the iterations when the significance level is 0.05, it can be done as
follows:

```r
summary(invertci.subsample2, alphas = 0.05)
#> Maximum number of iterations: 5
#> Tolerance level: 0.001
#>                                       
#> Significance level Confidence interval
#> 0.05                [0.33125, 0.66875]
#> 
#> Details:
#> 
#> <Confidence interval for significance level = 0.05>
#> 
#> === Iterations in constructing upper bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.60000            NA      0.60000   0.48000     FALSE
#> Right end pt.            NA       1.00000      1.00000   0.00000      TRUE
#> 1                   0.60000       1.00000      0.80000   0.00000      TRUE
#> 2                   0.60000       0.80000      0.70000   0.00000      TRUE
#> 3                   0.60000       0.70000      0.65000   0.15000     FALSE
#> 4                   0.65000       0.70000      0.67500   0.01000      TRUE
#> 5                   0.65000       0.67500      0.66250   0.05000     FALSE
#> Reason for termination: Reached maximum number of iterations
#> 
#> 
#> === Iterations in constructing lower bound:                                                                          
#> Iteration       Lower bound   Upper bound   Test point   p-value   Reject?
#> Left end pt.        0.00000            NA      0.00000   0.00000      TRUE
#> Right end pt.            NA       0.40000      0.40000   0.46000     FALSE
#> 1                   0.00000       0.40000      0.20000   0.00000      TRUE
#> 2                   0.20000       0.40000      0.30000   0.00000      TRUE
#> 3                   0.30000       0.40000      0.35000   0.05000     FALSE
#> 4                   0.30000       0.35000      0.32500   0.00000      TRUE
#> 5                   0.32500       0.35000      0.33750   0.00000      TRUE
#> Reason for termination: Reached maximum number of iterations
```

Parallel Programming
--------------------

The `lpinfer` package supports parallel programming through the `future`
package. The `future` package is freely available on CRAN and can be
installed directly via the `install.packages` command in R.

To implement parallel programming, users need to first install and load
the `future` package. Then, specify the number of workers to be used via
the `plan` function before running any functions from the `lpinfer`
package via the `plan` function after loading the `future` package.

Below is an example of running the `dkqs` function with parallel
programming by the `future` package:

```r
library(future)
plan(multisession, workers = 8)
set.seed(1)
dkqs.full1 <- dkqs(data = sampledata,
                   lpmodel = lpm.full,
                   beta.tgt = .375,
                   R = 1000,
                   tau = sqrt(log(N)/N),
                   solver = "gurobi",
                   progress = FALSE)
print(dkqs.full1)
```

The computational time for using multiple cores should be shorter than
using a single core for a large number of bootstraps. This is
illustrated by the example via the `dkqs` procedure below:

```r
library(future)
dkqs.args <- list(data = sampledata,
                  lpmodel = lpm.full,
                  beta.tgt = .375,
                  R = 5000,
                  tau = sqrt(log(N)/N),
                  solver = "gurobi",
                  progress = FALSE)

# Run dkqs with one worker
plan(multisession, workers = 1)
t1 <- Sys.time()
set.seed(1)
do.call(dkqs, dkqs.args)
#>  p-value: 0.2368
t2 <- Sys.time()
time1 <- t2 - t1

# Run dkqs with three workers
plan(multisession, workers = 3)
t1 <- Sys.time()
set.seed(1)
do.call(dkqs, dkqs.args)
#>  p-value: 0.2368
t2 <- Sys.time()
time3 <- t2 - t1

# Print the time used
print(sprintf("Time used with 1 core: %s", time1))
#> [1] "Time used with 1 core: 24.944885969162"
print(sprintf("Time used with 3 cores: %s", time3))
#> [1] "Time used with 3 cores: 12.7029490470886"
```

Note that there are different options available for the `plan` command.
The availability of the options may depend on the operating systems
used. Please refer to the instructions from the `future` pakage for more
details.

Help, Feature Requests and Bug Reports
--------------------------------------

Please post an issue on the [issues
page](https://github.com/conroylau/lpinfer/issues) of the GitHub
repository. We are happy to help.

References
----------

Deb, R., Y. Kitamura, J. K. H. Quah, and Stoye J. 2018. “Revealed Price
Preference: Theory and Empirical Analysis.” *Working Paper*.

Fang, Z., A. Santos, A. Shaikh, and A. Torgovitsky. 2020. *Working
Paper*.

Kamat, V. 2019. “Identification with Latent Choice Sets.” *Working
Paper*.

Manski, C. F. 1989. “Anatomy of the Selection Problem.” *The Journal of
Human Resources* 24: 343–60.
