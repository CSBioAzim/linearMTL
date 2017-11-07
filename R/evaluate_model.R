################################################################################
# Evaluate linear multi-task learning model.
################################################################################

#' Evaluate LinearMTL model.
#'
#' Compute training and test errors for the given indices as well as error
#' changes for excluding features. If task.grouping is supplied, identify top
#' regulators per group and plot clustering of the regression coefficients along
#' with the most important features. Save results to directory out.dir.
#'
#' @param X Column centered N by J input matrix of features common to all tasks.
#' @param task.specific.features Named list of features which are specific to
#'   each task. Each entry contains an N by J2 column-centered matrix for one
#'   particular task (where columns are features). List has to be ordered
#'   according to the columns of Y.
#' @param Y Column centered N by K output matrix for every task.
#' @param LMTL.model Linear multi-task learning model (list containing B and
#'   intercept).
#' @param train.idx Indices for the training set.
#' @param test.idx Indices for the test set.
#' @param task.grouping String vector of length K with group names for each
#'   task.
#' @param feature.names Feature names.
#' @param task.names Task names.
#' @param sd.threshold All features with error changes above sd.threshold times
#'   the standard deviation will be taken as important.
#' @param out.dir Output directory for results and plots.
#'
#' @export
EvaluateLinearMTModel <- function(X = NULL, task.specific.features = list(), Y, LMTL.model,
                                  train.idx, test.idx, task.grouping = NULL,
                                  feature.names = NULL, task.names = NULL,
                                  sd.threshold = 1.5, out.dir) {

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
  if (length(task.specific.features) > 0) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  if (nrow(LMTL.model$B) != J) {
    stop("B must have as many rows as there are features in X and task.specific.features!")
  }

  if (is.null(feature.names)) {
    feature.names <- 1:J
  }
  if (is.null(task.names)) {
    task.names <- 1:K
  }

  # compute predictions
  predictions <- MTPredict(LMTL.model = LMTL.model, X = X,
                           task.specific.features = task.specific.features)
  train.pred <- predictions[train.idx, ]
  test.pred <- predictions[test.idx, ]

  # compute training and test error
  err.out <- matrix(0, nrow = 2, ncol = 3)
  dimnames(err.out) <- list(c("train", "test"), c("mse", "cor", "r2"))

  train.residuals <- (Y[train.idx, ] - train.pred)^2
  test.residuals <- (Y[test.idx, ] - test.pred)^2

  train.mse <- colMeans(train.residuals)
  test.mse <- colMeans(test.residuals)

  train.cor.spearman <- diag(cor(train.pred, Y[train.idx, ], method = "spearman"))
  train.cor.pearson <- diag(cor(train.pred, Y[train.idx, ], method = "pearson"))

  test.cor.spearman <- diag(cor(test.pred, Y[test.idx, ], method = "spearman"))
  test.cor.pearson <- diag(cor(test.pred, Y[test.idx, ], method = "pearson"))

  train.r2 <- 1 - colSums(train.residuals) / colSums((Y[train.idx, ] - colMeans(Y[train.idx, ]))^2)
  test.r2 <- 1 - colSums(test.residuals) / colSums((Y[test.idx, ] - colMeans(Y[test.idx, ]))^2)

  err.out["train", "mse"] <- mean(train.mse)
  err.out["train", "cor"] <- mean(train.cor.spearman)
  err.out["train", "r2"] <- mean(train.r2)

  err.out["test", "mse"] <- mean(test.mse)
  err.out["test", "cor"] <- mean(test.cor.spearman)
  err.out["test", "r2"] <- mean(test.r2)

  # print matrix
  print(err.out)

  # save
  write.table(err.out, file.path(out.dir, "error.txt"),
              quote = FALSE, sep = '\t')


  task.by.task.statistics <- rbind(train.mse, test.mse,
                                   train.cor.pearson, test.cor.pearson,
                                   train.cor.spearman, test.cor.spearman)
  colnames(task.by.task.statistics) <- task.names
  saveRDS(task.by.task.statistics, file.path(out.dir, "tbt_stats.rds"))

  # compute error changes
  error.change <- matrix(0, J, K, dimnames = list(feature.names, task.names))

  print("Computing error changes ... ")
  # compute error change for excluding common features
  if (J1 > 0) {
    for (j in 1:J1) {
      feature.contribution <- X[train.idx, j, drop = FALSE] %*% LMTL.model$B[j, , drop = FALSE] + LMTL.model$intercept
      error.change[j, ] <- colSums(((train.pred - feature.contribution) - Y[train.idx, ])^2 - train.residuals)
    }
  }
  # compute error change for excluding task specific features
  if (length(task.specific.features) > 0) {
    for (j in 1:J2) {
      feature.contribution <- sapply(1:K, FUN = function(x){
        task.specific.features[[x]][train.idx, j, drop = FALSE] %*% LMTL.model$B[J1 + j, x, drop = FALSE] + LMTL.model$intercept[x]
      })
      error.change[J1 + j, ] <- colSums(((train.pred - feature.contribution) - Y[train.idx, ])^2 - train.residuals)
    }
  }

  error.change <- error.change / nrow(Y[train.idx,])
  saveRDS(error.change, file.path(out.dir, "error_changes.rds"))

  if (!is.null(task.grouping)) {
    dimnames(LMTL.model$B) <- list(feature.names, task.names)
    names(task.grouping) <- task.names

    # find top regulators for each group and plot
    print("Identifying top regulators ... ")
    tasks.by.group <- tapply(names(task.grouping), task.grouping, list)
    candidate.regulators <- list ()
    target.regulation <- list ()
    candidate.regulators$Common <- rownames(error.change)

    pdf(sprintf('%s/Agg_errors.pdf', out.dir))

    for (group in names(tasks.by.group)) {
      agg.errors <- rowSums(error.change[, tasks.by.group[[group]], drop = FALSE])
      cutoff <- mean(agg.errors) + sd.threshold * sd(agg.errors)
      agg.errors <- sort(agg.errors, decreasing = TRUE)
      candidates <- names(agg.errors)[agg.errors >= cutoff]
      candidate.regulators[[group]] <- candidates
      candidate.regulators$Common <- intersect(candidates, candidate.regulators$Common)
      target.regulation[[group]] <- apply(X = LMTL.model$B[,tasks.by.group[[group]], drop = FALSE], MARGIN = 1,
                                          FUN = function (x) {
                                            ifelse (length (which (x > 0)) > length (which (x<0)), 'Up', 'Down')
                                          })[candidates]
      df <- data.frame(x = 1:length(agg.errors), y = sort(agg.errors), Cutoff = 'Below', stringsAsFactors = FALSE)
      df$Cutoff[df$y >= cutoff] <- 'Above'
      plot (0, 0, xlim = c(1, length(agg.errors)), ylim = range(agg.errors), type = 'n',
            xlab = 'Regulators', ylab = 'Aggregate Error Changes', main = sprintf("%s", group))
      points(df$x[df$Cutoff == 'Below'], df$y[df$Cutoff == 'Below'], col = '#B2B2B2', pch = 20, cex = 1.5)
      points(df$x[df$Cutoff != 'Below'], df$y[df$Cutoff != 'Below'], pch = 20, cex = 1.5, col = '#EF264D')
      lines(c(-1000, 100000), c(cutoff, cutoff), lwd = 2)
      text(500, cutoff, 'Cutoff for candidate regulators\n(mean + 1.5*sd)', cex = 0.6, pos = 3)
    }

    dev.off ()

    if (sum(colSums(LMTL.model$B == 0) == J) == 0) {
      print("Clustering coefficient matrix ... ")
      pdf(sprintf("%s/Model_clustering.pdf", out.dir))

      # extract the selected features and overlay the dendrogram
      coefficient.dendrogram <- as.dendrogram(hclust(as.dist(1 - cor(LMTL.model$B, method = 'pearson'))))
      regulators <- unique(unlist(candidate.regulators))
      if (length(regulators) < 2) {
        regulators.to.plot <- union(feature.names[1:2], regulators)
      } else {
        regulators.to.plot <- regulators
      }
      PlotCustomHeatmap(matrix = as.matrix(LMTL.model$B[regulators.to.plot,]),
                        task.grouping = task.grouping,
                        Colv = coefficient.dendrogram)
      dev.off ()

    } else {
      print("No clustering of coefficient matrix due to constant regression vectors.")
    }

    # write candidate regulators and their target regulation to tab separated file
    max.length <- max(sapply(candidate.regulators, length))
    regulators.and.targets <- matrix ('',
                                      nrow = max.length,
                                      ncol = length(candidate.regulators) * 3)
    index <- 1
    for (group in names(tasks.by.group)) {
      num.of.regs <- length(candidate.regulators[[group]])
      if (num.of.regs > 0) {
        regulators.and.targets[1:num.of.regs, index]  <- candidate.regulators[[group]]
        regulators.and.targets[1:num.of.regs, index + 1]  <- target.regulation[[group]]
      } else {
        print(sprintf("Warning: Could not identify candidate regulators for group %s!", group))
        print("         Consider decreasing sd.threshold.")
      }
      index <- index + 3
    }
    if (length(candidate.regulators$Common) > 0) {
      regulators.and.targets[1:length(candidate.regulators$Common), index] <- candidate.regulators$Common
    }
    colnames(regulators.and.targets) <- rep('', ncol(regulators.and.targets))
    colnames(regulators.and.targets)[seq(1, ncol(regulators.and.targets), 3)] <- c(names(tasks.by.group), 'Common')

    # save
    write.table(regulators.and.targets, file.path(out.dir, "regulators_and_targets.txt"),
                quote = FALSE, sep = '\t')
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
#' @param Y Column centered N by K output matrix for every task.
#' @param LMTL.model.list List of LinearMTL models.
#' @param train.idx.by.cluster List of training indices per cluster.
#' @param test.idx.by.cluster List of test indices per cluster.
#' @param task.names Task names.
#' @param out.dir Output directory for results and plots.
#'
#' @export
EvaluateClusteredLinearMTModel <- function(X = NULL, task.specific.features = list(), Y, LMTL.model.list,
                                           train.idx.by.cluster, test.idx.by.cluster,
                                           task.names = NULL, out.dir) {

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
  if (length(task.specific.features) > 0) {
    if (nrow(task.specific.features[[1]]) != nrow(Y)) {
      stop("Task specific feature matrices and Y must have the same number of rows!")
    }
    J2 <- ncol(task.specific.features[[1]])
  }

  # input / output dimensions
  N <- nrow(Y)
  K <- ncol(Y)
  J <- J1 + J2

  if (is.null(task.names)) {
    task.names <- 1:K
  }

  predictions <- matrix(0, N, K)
  for (k in seq_along(LMTL.model.list)) {
    # compute predictions for every cluster
    train.idx <- train.idx.by.cluster[[k]]
    test.idx <- test.idx.by.cluster[[k]]
    LMTL.model <- LMTL.model.list[[k]]
    predictions[train.idx, ] <- MTPredict(LMTL.model = LMTL.model, X = X[train.idx, ],
                                          task.specific.features = lapply(task.specific.features,
                                                                          FUN = function(x){x[train.idx, ]}))
    predictions[test.idx, ] <- MTPredict(LMTL.model = LMTL.model, X = X[test.idx, ],
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
  err.out <- matrix(0, nrow = 2, ncol = 3)
  dimnames(err.out) <- list(c("train", "test"), c("mse", "cor", "r2"))

  train.residuals <- (Y[train.idx, ] - train.pred)^2
  test.residuals <- (Y[test.idx, ] - test.pred)^2

  train.mse <- colMeans(train.residuals)
  test.mse <- colMeans(test.residuals)

  train.cor.spearman <- diag(cor(train.pred, Y[train.idx, ], method = "spearman"))
  train.cor.pearson <- diag(cor(train.pred, Y[train.idx, ], method = "pearson"))

  test.cor.spearman <- diag(cor(test.pred, Y[test.idx, ], method = "spearman"))
  test.cor.pearson <- diag(cor(test.pred, Y[test.idx, ], method = "pearson"))

  train.r2 <- 1 - colSums(train.residuals) / colSums((Y[train.idx, ] - colMeans(Y[train.idx, ]))^2)
  test.r2 <- 1 - colSums(test.residuals) / colSums((Y[test.idx, ] - colMeans(Y[test.idx, ]))^2)

  err.out["train", "mse"] <- mean(train.mse)
  err.out["train", "cor"] <- mean(train.cor.spearman)
  err.out["train", "r2"] <- mean(train.r2)

  err.out["test", "mse"] <- mean(test.mse)
  err.out["test", "cor"] <- mean(test.cor.spearman)
  err.out["test", "r2"] <- mean(test.r2)

  # print / save
  print(err.out)

  # save
  write.table(err.out, file.path(out.dir, "error.txt"),
              quote = FALSE, sep = '\t')


  task.by.task.statistics <- rbind(train.mse, test.mse,
                                   train.cor.pearson, test.cor.pearson,
                                   train.cor.spearman, test.cor.spearman)
  colnames(task.by.task.statistics) <- task.names
  saveRDS(task.by.task.statistics, file.path(out.dir, "tbt_stats.rds"))
}




#' Plot customized heatmap.
#'
#' @param matrix Matrix to display.
#' @param scale See heatmap.2.
#' @param trace See heatmap.2
#' @param task.grouping Vector of strings with groups for each task.
#' @param dendro Dendrogram obtained from clustering the columns of the input matrix.
#' @param density.info See heatmap.2.
#' @param ... Other arguments passed to heatmap.2.
#'
#' @importFrom gplots heatmap.2
PlotCustomHeatmap <- function (matrix, scale = "none", trace = "none",
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
