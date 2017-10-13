################################################################################
# Evaluate linear multi-task learning model.
################################################################################

#' Evaluate LinearMTL model.
#'
#' Compute training and test errors for the given indices as well as error
#' changes for excluding features. If task.grouping is supplied, plot clustering
#' of the regression coefficients along with the most important features.
#'
#' @param X Column centered N by J input matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param B J by K matrix of regression coefficients, J = J1 + J2.
#' @param train.idx Indices for the training set.
#' @param test.idx Indices for the test set.
#' @param task.grouping String vector of length K with group names for each
#'   task.
#' @param feature.names Feature names.
#' @param task.names Task names.
#' @param sd.threshold All features with error changes above sd.threshold times
#'   the standard deviation will be taken as important.
#' @param out.dir Output directory.
#'
#' @export
EvaluateLinearMTModel <- function(X = NULL, task.specific.features = list(), Y, B,
                                  train.idx, test.idx, task.grouping = NULL,
                                  feature.names = NULL, task.names = NULL,
                                  sd.threshold = 1.5, out.dir = NULL) {

  # initialization and error checking
  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # check for shared features
  J1 <- 0
  if (!is.null(X)) {
    if (nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows!")
    }
    J1 <- ncol(X)
  }

  # check for task specific features
  J2 <- 0
  if (length(task.specific.features > 0)) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  if (nrow(B) != J) {
    stop("B must have as many rows as there are features in X and task.specific.features!")
  }

  # compute predictions
  predictions <- MTPredict(beta = B, X = X, task.specific.features = task.specific.features)
  train.pred <- predictions[train.idx, ]
  test.pred <- predictions[test.idx, ]

  # compute training and test error
  err.out <- matrix(0, nrow = 2, ncol = 2)
  dimnames(err.out) <- list(c("train", "test"), c("mse", "cor"))

  err.out["train", "mse"] <- MTComputeError(Y = Y[train.idx, ], beta = B, pred = train.pred)
  err.out["train", "cor"] <- MTComputeMeanCorrelation(Y = Y[train.idx, ], beta = B, pred = train.pred)

  err.out["test", "mse"] <- MTComputeError(Y = Y[test.idx, ], beta = B, pred = test.pred)
  err.out["test", "cor"] <- MTComputeMeanCorrelation(Y = Y[test.idx, ], beta = B, pred = test.pred)

  # print / save
  print(err.out)
  if (!is.null(out.dir)) {
    write.table(err.out, file = sprintf("%s/err.txt", out.dir),
                quote = FALSE, sep = '\t')
  }

  if (!is.null(out.dir) | !is.null(task.grouping)) {
    # compute error changes
    if (is.null(feature.names)) {
      feature.names <- 1:J
    }
    if (is.null(task.names)) {
      task.names <- 1:K
    }

    squared.diff <- MTComputeError(Y = Y[train.idx, ], beta = B, pred = train.pred, normalize = FALSE)
    error.change <- matrix(0, J, K, dimnames = list(feature.names, task.names))

    print("Computing error changes ... ")
    # compute error change for excluding common features
    for (j in 1:J1) {
      feature.contribution <- MTPredict(beta = B, X = X[train.idx, j, drop = FALSE])
      error.change[j, ] <- colSums(((train.pred - feature.contribution) - Y[train.idx, ])^2 - squared.diff)
    }
    # compute error change for excluding task specific features
    if (length(task.specific.features) > 0) {
      for (j in 1:J2) {
        feat.task.specific.features <- lapply(task.specific.features, FUN = function(x){x[train.idx, j, drop = FALSE]})
        feature.contribution <- MTPredict(beta = B, task.specific.features = feat.task.specific.features)
        error.change[J1 + j, ] <- colSums(((train.pred - feature.contribution) - Y[train.idx, ])^2 - squared.diff)
      }
    }

    error.change <- error.change / nrow(Y[train.idx,])
    # output error changes
    if (!is.null(out.dir)) {
      write.table(error.change,
                  file = sprintf("%s/Feature_selection.txt", out.dir),
                  quote = FALSE, sep = '\t')
    }
  }

  if (!is.null(task.grouping)) {
    dimnames(B) <- list(feature.names, task.names)
    names(task.grouping) <- task.names

    # find top regulators for each group and plot
    print("Identifying top regulators ... ")
    tasks.by.group <- tapply(names(task.grouping), task.grouping, list)
    candidate.regulators <- list ()
    target.regulation <- list ()
    candidate.regulators$Common <- rownames(error.change)

    if (!is.null(out.dir)) {
      pdf(sprintf('%s/Agg_errors.pdf', out.dir))
    }

    for (group in names(tasks.by.group)) {
      agg.errors <- rowSums(error.change[, tasks.by.group[[group]], drop = FALSE])
      cutoff <- mean(agg.errors) + sd.threshold * sd(agg.errors)
      agg.errors <- sort(agg.errors, decreasing = TRUE)
      candidates <- names(agg.errors)[agg.errors >= cutoff]
      candidate.regulators[[group]] <- candidates
      candidate.regulators$Common <- intersect(candidates, candidate.regulators$Common)
      target.regulation[[group]] <- apply(X = B[,tasks.by.group[[group]], drop = FALSE],MARGIN = 1,
                                          FUN = function (x) {ifelse (length (which (x > 0)) > length (which (x<0)), 'Up', 'Down') })[candidates]
      df <- data.frame(x = 1:length(agg.errors), y = sort(agg.errors), Cutoff = 'Below', stringsAsFactors = FALSE)
      df$Cutoff[df$y >= cutoff] <- 'Above'
      plot (0, 0, xlim = c(1, length(agg.errors)), ylim = range(agg.errors), type = 'n',
            xlab = 'Regulators', ylab = 'Aggregate Error Changes', main = sprintf("%s", group))
      points(df$x[df$Cutoff == 'Below'], df$y[df$Cutoff == 'Below'], col = '#B2B2B2', pch = 20, cex = 1.5)
      points(df$x[df$Cutoff != 'Below'], df$y[df$Cutoff != 'Below'], pch = 20, cex = 1.5, col = '#EF264D')
      lines(c(-1000, 100000), c(cutoff, cutoff), lwd = 2)
      text(500, cutoff, 'Cutoff for candidate regulators\n(mean + 1.5*sd)', cex = 0.6, pos = 3)
    }
    if (!is.null(out.dir)) {
      dev.off ()
    }

    print("Clustering coefficient matrix ... ")
    if (!is.null(out.dir)) {
      pdf(sprintf("%s/Model_clustering.pdf", out.dir))
    }

    # extract the selected features and overlay the dendrogram
    coefficient.dendrogram <- as.dendrogram(hclust(as.dist(1 - cor(B, method = 'pearson'))))
    regulators <- unique(unlist(candidate.regulators))
    if (length(regulators) < 2) {
      regulators.to.plot <- union(feature.names[1:2], regulators)
    } else {
      regulators.to.plot <- regulators
    }
    PlotCustomHeatmap(matrix = as.matrix(B[regulators.to.plot,]),
                      task.grouping = task.grouping,
                      Colv = coefficient.dendrogram)
    if (!is.null(out.dir)) {
      dev.off ()
    }

    # write candidate regulators and their target regulation to tab separated file
    max.length <- max(sapply(candidate.regulators, length))
    output <- matrix ('', nrow = max.length, ncol = length(candidate.regulators) * 3)
    index <- 1
    for (group in labels(tasks.by.group)) {
      output[1:length(candidate.regulators[[group]]), index]  <- candidate.regulators[[group]]
      output[1:length(candidate.regulators[[group]]), index + 1]  <- target.regulation[[group]]
      index <- index + 3
    }
    if (length(candidate.regulators$Common) > 0) {
      output[1:length(candidate.regulators$Common), index] <- candidate.regulators$Common
    }
    colnames(output) <- rep('', ncol(output))
    colnames(output)[seq(1, ncol(output), 3)] <- c(names(tasks.by.group), 'Common')

    if (!is.null(out.dir)) {
      write.table(output, file = sprintf("%s/Candidate_regulators.txt", out.dir),
                  quote = FALSE, sep = '\t', row.names = FALSE)
    }
  }
}

#' Evaluate clustered LinearMTL model.
#'
#' Compute training and test errors for the given indices.
#'
#' @param X Column centered N by J input matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param B.list List of J by K matrices of regression coefficients, J = J1 + J2.
#' @param train.idx.by.cluster List of training indices per cluster.
#' @param test.idx.by.cluster List of test indices per cluster.
#' @param out.dir Output directory.
#'
#' @export
EvaluateClusteredTreeGuidedGroupLasso <- function(X = NULL, task.specific.features = list(), Y, B.list,
                                                  train.idx.by.cluster, test.idx.by.cluster, out.dir = NULL) {

  # initialization and error checking
  if (is.null(X) & (length(task.specific.features) == 0)) {
    stop("No input data supplied.")
  }

  # check for shared features
  J1 <- 0
  if (!is.null(X)) {
    if (nrow(X) != nrow(Y)) {
      stop("X and Y must have the same number of rows!")
    }
    J1 <- ncol(X)
  }

  # check for task specific features
  J2 <- 0
  if (length(task.specific.features > 0)) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  predictions <- matrix(0, N, K)
  for (k in seq_along(B.list)) {
    # compute predictions for every cluster
    train.idx <- train.idx.by.cluster[[k]]
    test.idx <- test.idx.by.cluster[[k]]
    B <- B.list[[k]]

    predictions[train.idx, ] <- MTPredict(beta = B, X = X[train.idx, ],
                                          task.specific.features = lapply(task.specific.features,
                                                                          FUN = function(x){x[train.idx, ]}))
    predictions[test.idx, ] <- MTPredict(beta = B, X = X[test.idx, ],
                                         task.specific.features = lapply(task.specific.features,
                                                                         FUN = function(x){x[test.idx, ]}))
  }

  # combine all training sets
  train.idx <- unlist(train.idx.by.cluster)
  # combine all test sets
  test.idx <- unlist(test.idx.by.cluster)

  train.pred <- predictions[train.idx, ]
  test.pred <- predictions[test.idx, ]

  # compute training and test error
  err.out <- matrix(0, nrow = 2, ncol = 2)
  dimnames(err.out) <- list(c("train", "test"), c("mse", "cor"))

  err.out["train", "mse"] <- MTComputeError(Y = Y[train.idx, ], beta = B, pred = train.pred)
  err.out["train", "cor"] <- MTComputeMeanCorrelation(Y = Y[train.idx, ], beta = B, pred = train.pred)

  err.out["test", "mse"] <- MTComputeError(Y = Y[test.idx, ], beta = B, pred = test.pred)
  err.out["test", "cor"] <- MTComputeMeanCorrelation(Y = Y[test.idx, ], beta = B, pred = test.pred)

  # print / save
  print(err.out)
  if (!is.null(out.dir)) {
    write.table(err.out, file = sprintf("%s/err.txt", out.dir), quote = FALSE, sep = '\t')
  }
}




PlotCustomHeatmap <- function (matrix, scale = "none", trace = "none", RowV = TRUE,
                               task.grouping, dendro = "col", density.info = 'none', ...) {

  # heatmap colors
  col <- rep ("", 25)
  col[15:11] <- c('#e8e8e8', '#e0e0e0', '#d8d8d8', '#d0d0d0', '#c8c8c8')
  col[1:10] <- rev(c("#02E1F6","#04D1EC","#08B2DA","#0B96C7","#0D7CB4","#0E64A2","#0F4F8F","#0F3D7C","#0F2D69","#0E1F57"))
  col[16:25] <- c("#F4DE00","#ECCD07","#E3BD0D","#DBAD12","#CA901C","#C28221","#B97524","#B16928","#A85D2B","#A0522D")

  # assign a color to each group
  unique.types <- unique(task.grouping)
  colors.list <- rainbow(length(unique.types))
  names(colors.list) <- unique.types

  ColSideColors <- rep(colors.list[1], ncol(matrix))
  for (type in 1:length(unique.types)) {
    ColSideColors[task.grouping == unique.types[type]] <- colors.list[type]
  }

  ## Draw the heatmap
  heatmap.res <- gplots::heatmap.2(matrix, scale = scale, trace = trace, cexCol = 0.6, cexRow = 0.6,
                                   col = col, ColSideColors = ColSideColors, dendro = dendro,
                                   density.info = density.info, symbreaks = TRUE, ...)

  legend("bottomleft", as.vector(unique.types), cex = 0.7,
         fill = colors.list[1:length(unique.types)], bty = "n")
}
