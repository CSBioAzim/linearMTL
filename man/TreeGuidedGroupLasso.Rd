% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-guided_group_lasso.R
\name{TreeGuidedGroupLasso}
\alias{TreeGuidedGroupLasso}
\title{Fit a tree-guided group lasso model (Kim et al. 2010).}
\usage{
TreeGuidedGroupLasso(X = NULL, task.specific.features = list(), Y, groups,
  weights, lambda, max.iter = 10000, epsilon = 1e-05, mu = NULL,
  mu.adapt = 1, cached.mats = NULL, MSE.Lipschitz = NULL, init.B = NULL,
  verbose = 1, standardize = TRUE, fit.intercept = TRUE,
  row.weights = NULL)
}
\arguments{
\item{X}{N by J1 matrix of features common to all tasks.}

\item{task.specific.features}{List of features which are specific to each
task. Each entry contains an N by J2 matrix for one particular task (where
columns are features). List has to be ordered according to the columns of
Y.}

\item{Y}{N by K output matrix for every task.}

\item{groups}{Binary V by K matrix determining group membership: Task k in
group v iff groups[v,k] == 1.}

\item{weights}{V dimensional vector with group weights.}

\item{lambda}{Regularization parameter.}

\item{max.iter}{(Optional) Maximum number of iterations.}

\item{epsilon}{(Optional) Desired accuracy. If error change drops below
epsilon, the algorithm terminates.}

\item{mu}{(Optional) Determines accuracy of smooth approximation to the group
penalty term. If NULL, mu will be determined based on desired accuracy
epsilon. However, this may lead to numerical issues and it is thus
recommended to tune mu by hand.}

\item{mu.adapt}{(Optional) Multiply mu with a factor of mu.adapt every
iteration. Default is no adaptation (mu.adapt = 1).}

\item{cached.mats}{Precomputed matrices as produced by PrepareMatrices.}

\item{MSE.Lipschitz}{(Optional) Lipschitz constant for MSE. If missing, use
maximum Eigenvalue of XTX.}

\item{init.B}{(Optional) (J+1) by K matrix with initializations for the
regression coefficients and intercept.}

\item{verbose}{(Optional) Integer in {0,1,2}. verbose = 0: No output. verbose
= 1: Print summary at the end of the optimization. verbose = 2: Print
progress during optimization.}

\item{standardize}{(Optional) Default is TRUE. Standardize data (using R
function scale()). Coefficients will be returned on original scale.}

\item{fit.intercept}{(Optional) Default is TRUE. Include intercept.}

\item{row.weights}{(Optional) Use weighted MSE. When cached.mats is supplied,
it is assumed that the rows of X and Y were already weighted appropriately.}
}
\value{
List containing
\item{lambda}{Regularization parameter used.}
\item{weights}{Node weights used.}
\item{B}{Final estimate of the regression coefficients.}
\item{intercept}{Final estimate of the intercept terms.}
\item{obj}{Final objective value.}
\item{early.termination}{Boolean indicating whether the algorithm exceeded
max.iter iterations.}
}
\description{
Fit a tree-guided group lasso model via smoothed proximal gradient descent
(Kim et al. 2012). May be trained on shared or task specific feature
matrices.
}
\seealso{
\code{\link{RunGroupCrossvalidation}}
}
