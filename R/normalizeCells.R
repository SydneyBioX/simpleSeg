#' Normalizes and transforms cell data in preparation for clustering
#' (accepts dataframe, SingleCellExperiment and SpatialExperiment).
#'
#' @param cells A Dataframe of SingleCellExperiment or SpatialExperiment
#'              containing cells and features to be normalized/transformed
#' @param markers A list containing the names of cell markers which will be
#'                normalized and/or transformed.
#' @param assayIn If input is a SingleCellExperiment or SpatialExperiment with
#'                multiple assays, specify the assay to be normalized and/or
#'                transformed.
#' @param assayOut If input is a SingleCellExperiment or SpatialExperiment, the
#'                 new of the normalized data.
#' @param imageID If input is a SingleCellExperiment or SpatialExperiment, this
#'                is the name of the image ID variable in order to stratify.
#'                cells correctly
#' @param transformation The transformation/s to be performed, default is NULL,
#'                       accepted values: 'asinh' and 'sqrt'.
#' @param method The normalization method/s to be performed, default is NULL,
#'               accepted values: 'mean', 'minMax', 'trim99', 'PC1'.
#' @param cores The number or cores for parallel processing.
#'
#' @return returns a dataframe with individual cells as rows and features as
#'         columns.
#'
#' @examples
#'
#' library(cytomapper)
#' data("pancreasSCE")
#' cells.normalized <- normalizeCells(
#'   cells = pancreasSCE,
#'   markers = c("CD99", "PIN", "CD8a", "CDH"),
#'   assayIn = "counts",
#'   assayOut = "normCounts",
#'   imageID = "ImageNb",
#'   transformation = "asinh",
#'   method = "trim99"
#' )
#'
#' @export normalizeCells
#' @rdname normalizeCells
#' @importFrom SummarizedExperiment assay colData assay<-
#' @importFrom methods is
#' @importFrom stats resid lm prcomp quantile
normalizeCells <- function(cells,
                           markers = NULL,
                           assayIn = NULL,
                           assayOut = "norm",
                           imageID = "imageID",
                           transformation = NULL,
                           method = NULL,
                           cores = 1) {
  valid_transformations <- c("asinh", "sqrt")
  bad <- setdiff(transformation, valid_transformations)
  if (length(bad) > 0) {
    stop(
      "`", bad[1], "` is not a valid transformation. ",
      "Choose one of: `sqrt`, `asinh`"
    )
  }
  valid_methods <- c("mean", "minMax", "trim99", "PC1")
  bad <- setdiff(method, valid_methods)
  if (length(bad) > 0) {
    stop(
      "`", bad[1], "` is not a valid method. ",
      "Choose one of: `mean`, `minMax`, `trim99` or `PC1`"
    )
  }
  perImage <- intersect(method, c("mean", "minMax", "trim99"))

  sce <- NULL
  ## handling sce and se
  if (is(cells, "SingleCellExperiment") || is(cells, "SpatialExperiment")) {
    sce <- cells
    values <- if (is.null(assayIn)) assay(sce) else assay(sce, assayIn)
    values <- as.matrix(values)
    ids <- colData(sce)[[imageID]]
    if (is.null(ids) && length(perImage) > 0) {
      stop(
        "imageID '", imageID, "' is not a column of colData(cells); ",
        "it is needed to normalise each image for method ",
        paste0("`", perImage, "`", collapse = ", "), "."
      )
    }
    if (is.null(markers)) markers <- rownames(values)
    if (is.null(markers)) markers <- seq_len(nrow(values))
    missing <- setdiff(markers, rownames(values))
    X <- t(values[markers, , drop = FALSE])
  } else {
    if (!imageID %in% colnames(cells) && length(perImage) > 0) {
      stop(
        "imageID '", imageID, "' is not a column of cells; ",
        "it is needed to normalise each image for method ",
        paste0("`", perImage, "`", collapse = ", "), "."
      )
    }
    ids <- cells[[imageID]]
    if (is.null(markers)) {
      markers <- colnames(cells)[!colnames(cells) %in% imageID]
    }
    missing <- setdiff(markers, colnames(cells))
    X <- as.matrix(cells[, markers, drop = FALSE])
  }
  if (is.character(markers) && length(missing) > 0) {
    stop("Markers not found: ", paste(missing, collapse = ", "), ".")
  }

  for (t in transformation) {
    X <- switch(t,
      "asinh" = asinh(X),
      "sqrt" = {
        if (min(X) < 0) {
          stop(paste0(
            "Can't use `sqrt` normailastion on negative markers.",
            " Consider using `asinh`."
          ))
        }
        sqrt(X)
      }
    )
  }

  if (length(method) > 0) {
    groups <- if (!is.null(ids)) split(seq_len(nrow(X)), ids)
    for (m in method) {
      X <- switch(m,
        "mean" = .perImage(X, groups, meandiv),
        "minMax" = .perImage(X, groups, minMax),
        "trim99" = .perImage(X, groups, trim99),
        "PC1" = PC1(X)
      )
    }
  }

  if (!is.null(sce)) {
    values[markers, ] <- t(X)
    assay(sce, assayOut) <- values
    return(sce)
  }

  cells[, markers] <- X
  cells
}

## Apply a normalisation to the cells (rows of X) of each image separately.
.perImage <- function(X, groups, FUN) {
  for (g in groups) {
    X[g, ] <- FUN(X[g, , drop = FALSE])
  }
  X
}

meandiv <- function(x) {
  markerMeans <- apply(x, 2, mean, 0.001)
  markerMeans[markerMeans <= 0] <- 1
  sweep(x, 2, markerMeans, "/")
}

trim99 <- function(x) {
  apply(x, 2, function(x) {
    q <- quantile(x, 0.99)
    if (q <= 0) q <- 1
    pmin(x, q)
  })
}

minMax <- function(x) {
  apply(x, 2, function(x) {
    x <- pmax(x - min(x), 0)
    m <- max(x)
    if (m <= 0) m <- 1
    x / m
  })
}

PC1 <- function(X) {
  pca <- prcomp(X)

  PC1 <- pca$x[, "PC1"]

  for (i in seq_len(ncol(X))) {
    y <- X[, i]
    fit <- lm(y ~ PC1)
    int <- fit$coefficients["(Intercept)"]
    X[, i] <- pmax(resid(fit) + int, 0)
    q <- quantile(X[, i], 0.99)
    if (q <= 0) q <- 1
    X[, i] <- pmin(X[, i], q) / q
  }
  X
}
