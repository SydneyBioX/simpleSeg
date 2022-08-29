#' Normalizes and transforms cell data in preparation for clustering (accepts dataframes, singleCellExperiments and SpatialExperiments)
#'
#' @param cells A Dataframe of singlecellexperement or spatialexperiment containing cells and features to be normalized/transformed
#' @param markers A list containing the names of cell markers which will be normalized/transformed
#' @param assayIn If input is a SCE or SE with multiple assays, specify the assay to be normalized/transformed
#' @param assayOut If input is a SCE or SE, the new of the normalized data.
#' @param imageID If input is a SCE or SE, this is the name of the image ID variable in order to stratify cells correctly
#' @param transformation The transformation/s to be performed, default is NULL, accepted values: 'asinh' and 'sqrt'
#' @param method The normalization method/s to be performed, default is NULL, accepted values: 'mean', 'minMax', 'trim99', 'PC1'
#' @param cores The number or cores for parallel processing
#'
#' @return returns a dataframe with individual cells as rows and features as columns
#'
#' @examples
#' 
#' library(cytomapper)
#' data("pancreasSCE")
#' cells.normalized <- normalizeCells(cells = pancreasSCE, markers = c('CD99', 'PIN', 'CD8a', 'CDH'),
#'  assayIn = 'counts', assayOut = 'normCounts', imageID = 'ImageNb', transformation = 'asinh', method = 'trim99')
#'
#' @export normalizeCells
#' @rdname normalizeCells
#' @importFrom SummarizedExperiment assay colData assay<-
normalizeCells <- function(cells,
                           markers=NULL,
                           assayIn=NULL,
                           assayOut="norm",
                           imageID="imageID",
                           transformation=NULL,
                           method=NULL,
                           cores=1) {
    
    
    sce <- NULL
    ## handling sce and se
    if (is(cells, "SingleCellExperiment") | is(cells, "SpatialExperiment")) {
        sce <- cells
        if (is.null(assayIn)) {
            cells <- as.data.frame(t(assay(sce)))
        } else {
            cells <- as.data.frame(t(assay(sce,assayIn)))
        }
        cells[[imageID]] <- colData(sce)[[imageID]]
    }
    
    if (is.null(markers)) {
        markers <- colnames(cells)[!colnames(cells)%in%imageID]
    }
    
    PC1 <- function(cells, markers, imageID) {
      pca <- prcomp(cells[,markers])
      
      PC1 <- pca$x[,"PC1"]
      
      
      for (i in markers) {
        y <- cells[,i]
        fit <- lm(y ~ PC1)
        int <- fit$coefficients["(Intercept)"]
        cells[,i] <- pmax(resid(fit) + int, 0)
        q <- quantile(cells[,i], 0.99)
        if (q <= 0) q <- 1
        cells[,i] <- pmin(cells[,i], q) / q
      }
      return(cells)
      
    }

    if (!is.null(transformation)) {
      for (i in seq_along(transformation)) {
        cells[, markers] <- switch(transformation[i],
                                   "asinh"=asinh(cells[, markers]),
                                   "sqrt"=sqrt(cells[,markers]))
      }
    }
    
    if (!is.null(method)) {
      for (i in seq_along(method)) {
        cells <- switch(method[i],
                        "mean"=meandiv(cells, markers, imageID),
                        "minMax"=minMax(cells, markers, imageID),
                        "trim99"=trim99(cells, markers, imageID),
                        "PC1"=PC1(cells, markers, imageID))
      }
    }
    
    if (is(sce, "SingleCellExperiment") | is(cells, "SpatialExperiment")) {
      assay(sce, assayOut) <- t(cells[colnames(cells) != imageID])
      return(sce)
    }
    
    cells
}

meandiv <- function(cells, markers, imageID) {
  for (i in unique(cells[[imageID]])) {
    markerMeans <- apply(cells[cells[[imageID]] == i, markers], 
                         2, 
                         mean, 
                         0.001)
    markerMeans[markerMeans <= 0] <- 1
    cells[cells[[imageID]] == i, markers] <-
      sweep(cells[cells[[imageID]] == i, markers],
            2, 
            markerMeans,
            "/")
  }
  cells
}

trim99 <- function(cells, markers, imageID) {
  for (i in unique(cells[[imageID]])) {
    cells[cells[[imageID]] == i, markers] <- 
      apply(cells[cells[[imageID]] == i, markers], 2, function(x) {
      q <- quantile(x, 0.99)
      if (q<=0) q <- 1
      pmin(x, q)
    })
  }
  
  cells
}

minMax <- function(cells, markers, imageID) {
  for (i in unique(cells[[imageID]])) {
    cells[cells[[imageID]] == i, markers] <- 
      apply(cells[cells[[imageID]] == i, markers], 2, function(x) {
      x <- pmax(x - min(x), 0) 
      m <- max(x)
      if (m <= 0) m <- 1
      x / m
    })
  }
  
  cells
}
