#' Perform simple segmentation
#'
#' @param image An image
#' @param BPPARAM A BiocParallelParam object.
#'
#' @return A list of image masks
#'
#' @examples
#' 
#' 1+1
#'
#' @export
#' @rdname simpleSeg
#' @importFrom BiocParallel SerialParam bplapply

simpleSeg <- function(image,
                          nucleus_index = 1,
                          size_selection = 10,
                          smooth = 1,
                          tolerance = 0.01,
                          ext = 1){
    
    
    
    nuc <- image[,,nucleus_index]
    
    # Hotspot filtering. Intensities greater than the 99 percentile are changed to be exactly the 99th percentile intensity. This has the effect of removing outliers.
    nuc[nuc > quantile(nuc, 0.99,na.rm = TRUE)] <- quantile(nuc,
                                                            0.99,
                                                            na.rm = TRUE)
    nuc <- nuc/max(nuc,na.rm = TRUE)
    
    nuc1 <- nuc
    nuc1[is.na(nuc)] <- 0
    
    nuc1 <- gblur(sqrt(nuc1),
                  smooth)
    
    
    # Otsu thresholding. 
    nth <- otsu((nuc1),
                range = c(0,1)) # thresholding on the sqrt intensities works better. 
    nmask = nuc1 >nth # the threshold is squared to adjust of the sqrt previously.
    
    
    
    #### Watershed to segment nucleus.
    nuc1 <- nuc # make nuc a matrix
    
    nuc1[!nmask|is.na(nmask)] <- 0 #get the intensities values in the nmask, and set other values to 0 that are not part of nmask.
    
    
    # Size selection
    nMaskLabel <- bwlabel(nmask) 
    tnuc1 <- table(nMaskLabel)
    
    
    nuc1[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
    
    nmask[nMaskLabel%in%names(which(tnuc1<=size_selection))] <- 0 # sizes less than 10 are set to zero.
    
    
    
    kern = makeBrush(5, shape='disc')
    
    
    disk_blur <- filter2(sqrt(nuc1), kern)
    #disk_blur <- disk_blur/
    
    #creating a distance matrix
    #distNuc <- EBImage::distmap(nuc1)
    
    # water shed to segment
    #nuc1 <- (nuc1>0)*1
    #nuc1 <- EBImage::distmap(nuc1)
    nmask1 <- watershed(disk_blur * nmask,
                        tolerance = tolerance,
                        ext = ext)
    #nMaskLabel <- bwlabel(nmask1) 
    #nMaskLabel <- colorLabels(nmask1) 
    #display(nMaskLabel)
    kern = makeBrush(3, shape='disc')
    cell1 = dilate(nmask1, kern)
    disk1 = cell1-nmask1 >0
    disk1 <- watershed(disk1)
    if(whole_cell){
        return(cell1)
    }
    else{
        return(nmask1)
    }
    
    
}
simpleSegParalell <- function(image,
                              nucleus_index = 1,
                              size_selection = 10,
                              smooth = 1,
                              tolerance = 0.01,
                              ext = 1,
                              whole_cell = TRUE,
                              cores = 50){
    output <- BiocParallel::bplapply(image, simpleSeg, nucleus_index = nucleus_index, tolerance = tolerance,  size_selection = size_selection, smooth = smooth, whole_cell = whole_cell,
                                       BPPARAM  = BiocParallel::MulticoreParam(workers = cores))
}

cytSeg <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           tolerance = 0.01,
                   kernSize = 3
){
    
    
    kern = makeBrush(kernSize, shape='disc')
    
    cell = dilate(nmask, kern)
    
    
    disk = cell - nmask > 0
    
    
    
    longImage_disk <- data.frame(apply(asinh(image),3, as.vector), disk = as.vector(disk))
    #longImage_background <- data.frame(apply(asinh(background),3, as.vector), disk = FALSE)
    #long_Image_all_disk <- longImage_disk %>% dplyr::filter(disk == TRUE)
    #long_Image_all_disk <- as.data.frame(long_Image_all_disk)
    
    
    #long_image_2 <- rbind(long_Image_all_disk, longImage_background)
    long_image_2 <- longImage_disk
    
    fit <- lm(disk ~ .-disk, data = long_image_2)
    
    cytpred <- nmask
    cytpred[] <- predict(fit, longImage_disk)
    cytpred <- cytpred - min(cytpred)
    cytpred <- cytpred/max(cytpred)
    
    cellTh <- otsu(cytpred,range = c(0,1))
    cell <- cytpred > cellTh
    
    cell <- cell + nmask > 0
    
    nuc_label <- bwlabel(nmask) 
    tnuc <- table(nuc_label)
    nmask[nuc_label%in%names(which(tnuc<=size_selection))] <- 0
    
    
    cmask4 <- propagate(cytpred, nmask, cell)
    justdisk <- propagate(disk, nmask, cell)
    
    #output<-list(cmask4,cmaskdisk, cell1, disk, justdisk)
    return(cmask4)
    
}
cytSegParalell <- function(nmask,
                           image,
                           size_selection = 5,
                           smooth = 1,
                           tolerance = 0.01,
                           kernSize = 3,
                           cores = 50){
    test.masks.cyt <- mcmapply(cytSeg, nmask, image, size_selection, mc.cores = 40)
}


normalize.cells <- function(cells,
                            markers,
                            transformation = NULL,
                            method = NULL){
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