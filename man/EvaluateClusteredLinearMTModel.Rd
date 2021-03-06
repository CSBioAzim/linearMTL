% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluate_model.R
\name{EvaluateClusteredLinearMTModel}
\alias{EvaluateClusteredLinearMTModel}
\title{Evaluate clustered LinearMTL model.}
\usage{
EvaluateClusteredLinearMTModel(X = NULL, task.specific.features = list(), Y,
  LMTL.model.list, train.idx.by.cluster, test.idx.by.cluster,
  task.names = NULL, out.dir)
}
\arguments{
\item{X}{Column centered N by J input matrix of features common to all tasks.}

\item{task.specific.features}{Named list of features which are specific to
each task. Each entry contains an N by J2 column-centered matrix for one
particular task (where columns are features). List has to be ordered
according to the columns of Y.}

\item{Y}{Column centered N by K output matrix for every task.}

\item{LMTL.model.list}{List of LinearMTL models.}

\item{train.idx.by.cluster}{List of training indices per cluster.}

\item{test.idx.by.cluster}{List of test indices per cluster.}

\item{task.names}{Task names.}

\item{out.dir}{Output directory for results and plots.}
}
\description{
Compute training and test errors for the given indices.
}
