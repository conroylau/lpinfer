context("Unit test for writing linear program in standard form")
rm(list = ls())

# ---------------- #
# Create the constraint matrices and rhs vectors
# ---------------- #
# Obs
Aobs <- matrix(c(1,1,2,3,3,7), nrow = 2, byrow = TRUE)
bobs <- matrix(c(6,20), nrow = 2, byrow = TRUE)
sobs <- matrix(rep("=", 2), nrow = 2, byrow = TRUE)

# Shp
Ashp <- matrix(c(10,0,0,2,1,2), nrow = 2, byrow = TRUE)
bshp <- matrix(c(1,10), nrow = 2, byrow = TRUE)
sshp <- matrix(c(">=", "<="), nrow = 2, byrow = TRUE)

# Tgt
Atgt <- matrix(c(1,2,3), nrow = 1, byrow = TRUE)
btgt <- c(9)
stgt <- c("=")

# lb and ub
xlb <- c(0.01,0.02,0.03)
xub <- c(1000,2000,3000)

# Get lpmodel.natural
lpmn <- lpmodel.natural(A.obs = Aobs,
                        A.shp = Ashp,
                        A.tgt = Atgt,
                        beta.obs = bobs,
                        beta.shp = bshp,
                        sense.shp = sshp,
                        x.lb = xlb,
                        x.ub = xub)
lpm <- standard.lpmodel(lpmn)

# ---------------- #
# Solve the model by Gurobi (not in standard form)
# ---------------- #
model1 <- list()

# Objective function
model1$obj <- c(1,2,3)

# Linear constraints
model1$A <- rbind(Aobs, Ashp, Atgt)
model1$rhs <- c(bobs, bshp, btgt)

# Model sense, lower bound and upper bound
model1$sense <- c(sobs, sshp, stgt)
model1$modelsense <- "min"
model1$lb <- xlb
model1$ub <- xub

# Get solution
params1 <- list(OutputFlag = 0, FeasibilityTol = 1e-9)
solution1 <- gurobi::gurobi(model1, params1)

# ---------------- #
# Solve the model by Gurobi (in standard form)
# ---------------- #
model2 <- list()

# Objective function (Extended)
model2$obj <- c(1,2,3, rep(0, ncol(lpm$A.obs) - 3))

# Linear constraints
model2$A <- rbind(lpm$A.obs, lpm$A.shp, lpm$A.tgt)
model2$rhs <- c(lpm$beta.obs, lpm$beta.shp, btgt)

# Model sense and lower bound
model2$modelsense <- "min"
model2$lb <- c(xlb, rep(0, ncol(lpm$A.obs) - 3))
model2$sense <- rep("=", nrow(model2$A))

# Get solution
params2 <- list(OutputFlag = 0, FeasibilityTol = 1e-9)
solution2 <- gurobi::gurobi(model2, params2)

# ---------------- #
# Run tests
# ---------------- #
# Objective value
test_that("Objective value",{
  expect_equal(solution1$objval, solution2$objval)
})

# Optimal point
test_that("Optimal point",{
  expect_equal(solution1$x, solution2$x[1:ncol(Ashp)])
})
