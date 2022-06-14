normalize.cells <- function(cells,
                            markers,
                            isSCE = FALSE,
                            assayName = NULL,
                            imageNb = NULL, #requires the image numbering varlable from sce list data to stratify the cells on
                            transformation = NULL,
                            method = NULL){
    
    #handeling sce and se
    if (isSCE == TRUE){
        if (is.null(assayName)){
            cellsdf <- data.frame(t(assay(cells)))
        }
        else{
            cellsdf <- data.frame(t(cells@assays@data@listData[[assayName]]))
        }
        cellsdf$imageID <- cells@colData@listData[[imageNb]]
        cells <- cellsdf
    }
    
    #Transformations
    if (transformation == "asinh"){
        cells[,markers] <- asinh(cells[,markers])
    }
    if (transformation == "sqrt"){
        cells[,markers] <- sqrt(cells[,markers])
    }
    if (transformation == "log"){
        cells[,markers] <- log10(cells[,markers])
    }
    #Methods
    if (method == "meandiv"){
        for(i in 1:length(unique(cells$imageID))){
            cells[cells$imageID==i,markers] <- sweep(cells[cells$imageID==i,markers], 2, apply(cells[cells$imageID==i,markers],2,mean, 0.2), "/")
        }
    }
    if (method == "99perc"){
        cells[,markers] <- data.frame(apply(cells[,markers],
                                            2, 
                                            function(x){
                                                q <- quantile(x,0.99)
                                                pmin(x,q)/q
                                            }))
    }
    if (method == "1stPC"){
        pca <- prcomp(cells[,markers])
        
        PC1 <- pca$x[,"PC1"]
        
        
        for(i in markers){
            y <- cells[,i]
            fit <- lm(y~PC1)
            int <- fit$coefficients["(Intercept)"]
            cells[,i] <- pmax(resid(fit) + int,0)
            q <- quantile(cells[,i],0.99)
            cells[,i] <- pmin(cells[,i],q)/q
        }
        
        
    }
    if (method == "scMerge"){
        #sc merge
        ncores <- 64
        use_bpparam <- BiocParallel::MulticoreParam(workers = ncores)
        use_bsparam <- BiocSingular::RandomParam()
        use_bnparam <- BiocNeighbors::AnnoyParam()
        # dat_sub <- dat[sample(nrow(dat), 500000), ]
        # ctl_genes <- rownames(sce)
        exprsMat <- t(cells[, markers])
        colnames(exprsMat) <- seq_len(ncol(exprsMat))
        scMerge_res <- scMerge2(exprsMat = exprsMat, #the exprs matrix to be normalised
                                batch = cells$imageID, # batch labels
                                cellTypes = NULL, # set NULL clustering will be performed within scMerge2... can also try the published cell type labels, which will match between the cell types
                                use_bpparam = use_bpparam,
                                use_bsparam = use_bsparam,
                                use_bnparam = use_bnparam,
                                ruvK = 2, # Number of unwanted variation to be removed
                                ctl = markers, # negative control genes
                                k_psuedoBulk = 5, # Number of pseudo bulk to be created for each cell type each batch
                                k_celltype = 20, # Number of neighbours when perform graph clustering
                                pseudoBulk_fn = create_pseudoBulk, # ways of constructing pseudo bulk
                                ncores = ncores,
                                chosen.hvg = markers, #Highly variable genes to be used to identify pseudo-replicates... since IMC has very few features, using all features.
                                cosineNorm = F,
                                return_subset = FALSE,
                                normalised = T)
        dat_norm <- as.data.frame(scMerge_res$newY)
        cells[, markers] <- dat_norm
    }
    return(cells)
}
