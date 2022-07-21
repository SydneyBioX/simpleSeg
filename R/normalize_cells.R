#' Normalizes and transforms cell data in preparation for clustering (accepts dataframes, singleCellExperiments and SpatialExperiments)
#'
#' @param cells A Dataframe of singlecellexperement or spatialexperiment containing cells and features to be normalized/transformed
#' @param markers A list containing the names of cell markers which will be normalized/transformed
#' @param assayIn If input is a SCE or SE with multiple assays, specify the assay to be normalized/transformed
#' @param assayOut If input is a SCE or SE, the new of the normalized data.
#' @param imageID If input is a SCE or SE, this is the name of the image ID variable in order to stratify cells correctly
#' @param transformation The transformation/s to be performed, default is NULL, accepted values: 'asinh', 'sqrt', 'log'
#' @param method The normalization method/s to be performed, default is NULL, accepted values: 'meandiv', '99perc', '1stPC'
#'
#' @return returns a dataframe with individual cells as rows and features as columns

#' @examples
#'
#' cells.normalized <- normalizeCells(cells = cells.sce, markers = c('SMA', 'CD44', 'CD45', 'cyt-19') isSCE = TRUE, assayName = 'cells', imageNb = 'ImageID', transformation = 'asinh', method = '99perc')
#'
#' @export normalizeCells
#' @rdname normalizeCells
#'
normalizeCells <- function(cells,
                           markers = NULL,
                           assayIn = NULL,
                           assayOut = "norm",
                           imageID = "imageID",
                           transformation = NULL,
                           method = NULL,
                           cores = 1) {
    
    
    sce <- NULL
    # handeling sce and se
    if (is(cells, "SingleCellExperiment")|is(cells, "SpatialExperiment")) {
        sce <- cells
        if (is.null(assayIn)) {
            cells <- as.data.frame(t(assay(sce)))
        } else {
            cells <- as.data.frame(t(assay(sce,assayIn)))
        }
        cells[[imageID]] <- colData(sce)[[imageID]]
    }
    
    if(is.null(markers)) {
        markers <- colnames(cells)[!colnames(cells)%in%imageID]
    }
    
    # Transformations
    #if (transformation == "asinh") {
    #    cells[, markers] <- asinh(cells[, markers])
    #}
    #if (transformation == "sqrt") {
    #    cells[, markers] <- sqrt(cells[, markers])
    #}
    #if (transformation == "log") {
    #    cells[, markers] <- log10(cells[, markers])
    #}
    
    # Methods
    # if ("meandiv" %in% method) {
    #     for (i in unique(cells$imageID)) {
    #         markerMeans <- apply(cells[cells[[imageID]] == i, markers], 2, mean, 0.001)
    #         markerMeans[markerMeans <= 0] <- 1
    #         cells[cells[[imageID]] == i, markers] <- sweep(cells[cells[[imageID]] == i, markers], 2, markerMeans,
    #                                                        "/")
    #     }
    # }
    # if ("99perc" %in% method) {
    #     for (i in unique(cells$imageID)) {
    #         cells[cells[[imageID]] == i, markers] <- apply(cells[cells[[imageID]] == i, markers], 2, function(x) {
    #             q <- quantile(x, 0.99)
    #             if(q<=0) q <- 1
    #             pmin(x, q) / q
    #         })
    #     }
    # }
    # 
    # if ("1stPC" %in% method) {
    #     pca <- prcomp(cells[, markers])
    #     
    #     PC1 <- pca$x[, "PC1"]
    #     
    #     
    #     for (i in markers) {
    #         y <- cells[, i]
    #         fit <- lm(y ~ PC1)
    #         int <- fit$coefficients["(Intercept)"]
    #         cells[, i] <- pmax(resid(fit) +
    #                                int, 0)
    #         q <- quantile(cells[, i], 0.99)
    #         if(q <= 0) q <- 1
    #         cells[, i] <- pmin(cells[, i], q) / q
    #     }
    #     
    #     
    # }
    # 
    # if ('scMerge' %in% method){
    #     if(!requireNamespace("scMerge", quietly = TRUE))
    #         stop("The package 'scMerge' could not be found. Please install it.")
    #     if(!requireNamespace("scMerge", quietly = TRUE))
    #         stop("The package 'scMerge' could not be found. Please install it.")
    #     if(!requireNamespace("scMerge", quietly = TRUE))
    #         stop("The package 'scMerge' could not be found. Please install it.")
    #     use_bpparam <- generateBPParam(cores)
    #     use_bsparam <- BiocSingular::RandomParam()
    #     use_bnparam <- BiocNeighbors::AnnoyParam()
    #     dat_sub <- dat[sample(nrow(dat), 500000),]
    #     ctl_genes <- rownames(sce)
    #     exprsMat <- t(cells[, markers])
    #     colnames(exprsMat) <- seq_len(ncol(exprsMat))
    #     scMerge_res <-
    #         scMerge2(
    #             exprsMat = exprsMat,
    #             #the exprs matrix to be normalised
    #             batch = cells[[imageID]],
    #             # batch labels
    #             cellTypes = NULL,
    #             # set NULL clustering will be performed within scMerge2... can also try the published cell type labels, which will matchetween the cell types
    #             use_bpparam = use_bpparam,
    #             use_bsparam = use_bsparam,
    #             use_bnparam = use_bnparam,
    #             ruvK = 2,
    #             # Number of unwanted variation to be removed
    #             ctl = markers,
    #             # negative control genes
    #             k_psuedoBulk = 5,
    #             # Number of pseudo bulk to be created for each cell type each batch
    #             k_celltype = 20,
    #             # Number of neighbours when performgraph clustering
    #             pseudoBulk_fn = create_pseudoBulk,
    #             # ways of onstructing pseudo bulk
    #             ncores = ncores,
    #             chosen.hvg = markers,
    #             #Highlyvariable genes to be used to identify pseudo-replicates... since IMC has  very few features, using all features.
    #             cosineNorm = F,
    #             return_subset = FALSE,
    #             normalised = T
    #         )
    #     dat_norm <- as.data.frame(scMerge_res$newY)
    #     cells[, markers] <- dat_norm
    # }
    
    
    #############################################FUNCTIONS########################################################
    
    ########################################transformations#####################################################
    # if (transformation == "asinh") {
    #   cells[, markers] <- asinh(cells[, markers])
    # }
    # if (transformation == "sqrt") {
    #   cells[, markers] <- sqrt(cells[, markers])
    # }
    # if (transformation == "log") {
    #   cells[, markers] <- log10(cells[, markers])
    # }
    # 
    ###################################################methods################################################
    meandiv <- function(cells, markers, imageID){
      for (i in unique(cells[[imageID]])) {
        markerMeans <- apply(cells[cells[[imageID]] == i, markers], 2, mean, 0.001)
        markerMeans[markerMeans <= 0] <- 1
        cells[cells[[imageID]] == i, markers] <- sweep(cells[cells[[imageID]] == i, markers], 2, markerMeans,
                                                       "/")
      }
      return(cells)
    }
    perc99 <- function(cells, markers, imageID) {
      for (i in unique(cells[[imageID]])) {
        cells[cells[[imageID]] == i, markers] <- apply(cells[cells[[imageID]] == i, markers], 2, function(x) {
          q <- quantile(x, 0.99)
          if(q<=0) q <- 1
          pmin(x, q) / q
        })
      }
      return(cells)
    }
    
    PC1 <- function(cells,markers, imageID) {
      pca <- prcomp(cells[, markers])
      
      PC1 <- pca$x[, "PC1"]
      
      
      for (i in markers) {
        y <- cells[, i]
        fit <- lm(y ~ PC1)
        int <- fit$coefficients["(Intercept)"]
        cells[, i] <- pmax(resid(fit) +
                             int, 0)
        q <- quantile(cells[, i], 0.99)
        if(q <= 0) q <- 1
        cells[, i] <- pmin(cells[, i], q) / q
      }
      return(cells)
      
    }
    
    Mergesc<- function(cells,markers, imageID){
      if(!requireNamespace("scMerge", quietly = TRUE))
        stop("The package 'scMerge' could not be found. Please install it.")
      if(!requireNamespace("scMerge", quietly = TRUE))
        stop("The package 'scMerge' could not be found. Please install it.")
      if(!requireNamespace("scMerge", quietly = TRUE))
        stop("The package 'scMerge' could not be found. Please install it.")
      use_bpparam <- generateBPParam(cores)
      use_bsparam <- BiocSingular::RandomParam()
      use_bnparam <- BiocNeighbors::AnnoyParam()
      dat_sub <- dat[sample(nrow(dat), 500000),]
      ctl_genes <- rownames(sce)
      exprsMat <- t(cells[, markers])
      colnames(exprsMat) <- seq_len(ncol(exprsMat))
      scMerge_res <-
        scMerge2(
          exprsMat = exprsMat,
          #the exprs matrix to be normalised
          batch = cells[[imageID]],
          # batch labels
          cellTypes = NULL,
          # set NULL clustering will be performed within scMerge2... can also try the published cell type labels, which will matchetween the cell types
          use_bpparam = use_bpparam,
          use_bsparam = use_bsparam,
          use_bnparam = use_bnparam,
          ruvK = 2,
          # Number of unwanted variation to be removed
          ctl = markers,
          # negative control genes
          k_psuedoBulk = 5,
          # Number of pseudo bulk to be created for each cell type each batch
          k_celltype = 20,
          # Number of neighbours when performgraph clustering
          pseudoBulk_fn = create_pseudoBulk,
          # ways of onstructing pseudo bulk
          ncores = ncores,
          chosen.hvg = markers,
          #Highlyvariable genes to be used to identify pseudo-replicates... since IMC has  very few features, using all features.
          cosineNorm = F,
          return_subset = FALSE,
          normalised = T
        )
      dat_norm <- as.data.frame(scMerge_res$newY)
      cells[, markers] <- dat_norm
      
      return(cells)
    }
    if(is.null(transformation) == FALSE){
      for (i in 1:length(transformation)){
        cells[, markers] <- switch(i,
                        "asinh" = asinh(cells[, markers]),
                        "sqrt" = sqrt(cells[,markers]))
                        #"log" = log10(cells[,markers]))
      }
    }
    
    if (is.null(method) == FALSE){
      for (i in 1:length(method)){
        cells <- switch(i,
                        "meandiv" = meandiv(cells, markers, imageID),
                        "perc99" = perc99(cells, markers, imageID),
                        "PC1" = PC1(cells, markers, imageID),
                        "Mergesc" = Mergesc(cells, markers, imageID)
        )
      }
    }
    
    if(is(sce, "SingleCellExperiment")|is(cells, "SpatialExperiment")){
      assay(sce, assayOut) <- t(cells[colnames(cells)!=imageID])
      return(sce)
    }
    
    
    return(cells)
}
